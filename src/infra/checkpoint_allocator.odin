// SPDX-FileCopyrightText: 2025 Rowan Apps
// SPDX-License-Identifier: MIT
/*
 Stack allocator based on checkpoints.
 -------------------------------------
 USE CASES:
    - Memory is allocated and freed in a LIFO pattern.
    - Need to free multiple allocations at once, but not nessecarily all previous allocations.
    - Not expected to use the same allocator across threads.
    - No need for error handling (other than allocator errors to match the interface).

Uses `assert` and relies on Odin's bounds checking to catch logic bugs with check & rewind.
Does implement Odin allocator interface so can be used with `make()`, `new()` etc.
*/
package infra

import sa "core:container/small_array"
import "core:mem"
import "core:mem/virtual"
import "core:slice"

MAX_CHECKPOINTS :: 32

Checkpoint_Allocator :: struct {
	backing:     virtual.Arena,
	data:        []byte,
	offset:      int,
	checkpoints: sa.Small_Array(MAX_CHECKPOINTS, int),
}

Checkpoint :: int

// Create a new checkpoint allocator, uses a virtual arena under the hood for backing.
ca_init :: proc(ca: ^Checkpoint_Allocator) -> mem.Allocator_Error {
	virtual.arena_init_growing(&ca.backing) or_return
	ca.data = make([]byte, mem.DEFAULT_PAGE_SIZE * 4, virtual.arena_allocator(&ca.backing)) or_return
	return nil
}

// Create a new checkpoint allocator backed with provided buffer.
ca_init_buffered :: proc(ca: ^Checkpoint_Allocator, buffer: []byte) -> mem.Allocator_Error {
	virtual.arena_init_buffer(&ca.backing, buffer) or_return
	ca.data = make([]byte, mem.DEFAULT_PAGE_SIZE * 4, virtual.arena_allocator(&ca.backing)) or_return
	return nil
}

// Free memory back to OS, if initialized with a buffer does not delete provided buffer.
ca_deinit :: proc(ca: ^Checkpoint_Allocator) {
	virtual.arena_destroy(&ca.backing)
}

// New checkpoint
ca_check :: proc(ca: ^Checkpoint_Allocator) -> (c: Checkpoint) {
	assert(sa.push_back(&ca.checkpoints, ca.offset), "Failed to create checkpoint, too many created.")
	return sa.len(ca.checkpoints) - 1
}

// Rewind and remove most recent checkpoint
ca_rewind :: proc(ca: ^Checkpoint_Allocator) {
	checkpoint := sa.pop_back(&ca.checkpoints)
	ca.offset = checkpoint
}

// Rewind and keep most recent checkpoint
ca_rewind_keep :: proc(ca: ^Checkpoint_Allocator) {
	checkpoint := ca.checkpoints.data[ca.checkpoints.len - 1]
	ca.offset = checkpoint
}

// Rewind to specific checkpoint and discard it
ca_rewind_to :: proc(ca: ^Checkpoint_Allocator, checkpoint: Checkpoint) {
	ca.offset = ca.checkpoints.data[checkpoint]
	ca.checkpoints.len = checkpoint // removes the checkpoint
}

// Rewind to a specific checkpoint and keep it
ca_rewind_to_keep :: proc(ca: ^Checkpoint_Allocator, checkpoint: Checkpoint) {
	ca.offset = ca.checkpoints.data[checkpoint]
	ca.checkpoints.len = checkpoint + 1 // + 1 to keep the current checkpoint
}

// Resize the data buffer, only call manually if you want to provide more size up front, cannot shrink the backing data.
ca_resize_backing :: proc(ca: ^Checkpoint_Allocator, new_size: int) -> mem.Allocator_Error {
	if new_size <= len(ca.data) { return nil }

	// Since we're using virtual arena, resize should preserve pointers
	backing_allocator := virtual.arena_allocator(&ca.backing)
	new_data, err := mem.resize(raw_data(ca.data), slice.size(ca.data), new_size, allocator = backing_allocator)
	if err != nil { return err }

	ca.data = slice.from_ptr(cast(^u8)new_data, new_size)
	return nil
}

// Create an allocator interface from the checkpoint allocator
ca_allocator :: proc(ca: ^Checkpoint_Allocator) -> mem.Allocator {
	return mem.Allocator{procedure = ca_allocator_proc, data = ca}
}

// Odin allocator interface implementation
ca_allocator_proc :: proc(
	allocator_data: rawptr,
	mode: mem.Allocator_Mode,
	size, alignment: int,
	old_memory: rawptr,
	old_size: int,
	loc := #caller_location,
) -> (
	[]byte,
	mem.Allocator_Error,
) {
	ca := cast(^Checkpoint_Allocator)allocator_data

	switch mode {
	case .Alloc, .Alloc_Non_Zeroed:
		aligned_offset := mem.align_forward_int(ca.offset, alignment)
		end_offset := aligned_offset + size

		// Check if we need to resize backing
		if end_offset > len(ca.data) {
			new_size := max(len(ca.data) * 2, end_offset) // TODO: doubling might be a bit agressive, maybe align end offset to a page boundary and allocate that.
			if err := ca_resize_backing(ca, new_size); err != nil { return nil, err }
		}

		result := ca.data[aligned_offset:end_offset]
		ca.offset = end_offset

		if mode == .Alloc { mem.zero_slice(result) }

		return result, nil
	case .Free:
		return nil, nil // no op on free.
	case .Free_All:
		ca.offset = 0
		sa.clear(&ca.checkpoints)
		return nil, nil
	case .Resize, .Resize_Non_Zeroed:
		if old_memory == nil { return ca_allocator_proc(allocator_data, .Alloc if mode == .Resize else .Alloc_Non_Zeroed, size, alignment, nil, 0, loc) }

		old_ptr := cast(uintptr)old_memory
		data_start := cast(uintptr)raw_data(ca.data)

		if old_ptr + uintptr(old_size) == data_start + uintptr(ca.offset) {
			aligned_offset := mem.align_forward_int(ca.offset - old_size, alignment)
			new_end_offset := aligned_offset + size

			if new_end_offset > len(ca.data) {
				new_backing_size := max(len(ca.data) * 2, new_end_offset)
				if err := ca_resize_backing(ca, new_backing_size); err != nil { return nil, err }
			}

			ca.offset = new_end_offset
			result := ca.data[aligned_offset:new_end_offset]

			if mode == .Resize && size > old_size { mem.zero_slice(result[old_size:]) }

			return result, nil
		} else {
			new_memory, err := ca_allocator_proc(allocator_data, .Alloc if mode == .Resize else .Alloc_Non_Zeroed, size, alignment, nil, 0, loc)
			if err != nil { return nil, err }

			copy_size := min(old_size, size)
			copy(new_memory[:copy_size], slice.from_ptr(cast(^u8)old_memory, old_size)[:copy_size])

			return new_memory, nil
		}
	case .Query_Features:
		set := (^mem.Allocator_Mode_Set)(old_memory)
		if set != nil {
			set^ = {.Alloc, .Alloc_Non_Zeroed, .Free, .Free_All, .Query_Features, .Resize, .Resize_Non_Zeroed}
		}
		return nil, nil
	case .Query_Info:
		return nil, .Mode_Not_Implemented
	}
	return nil, .Mode_Not_Implemented
}