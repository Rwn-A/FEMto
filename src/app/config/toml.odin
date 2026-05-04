package config

import "base:intrinsics"
import "core:fmt"
import "core:log"
import "core:mem"
import "core:reflect"
import "core:strings"


import "../../../vendor/toml"


@(private)
parse_tag :: proc(tag: string) -> (key: string, required: bool) {
	comma := strings.index_byte(tag, ',')
	if comma == -1 {return "", tag == "required"}
	return tag[comma + 1:], tag[:comma] == "required"
}

@(private)
assign_int :: proc(dest: any, i: i64) {
	v := reflect.any_core(dest)
	switch &dst in v {
	case i8:
		dst = i8(i)
	case i16:
		dst = i16(i)
	case i32:
		dst = i32(i)
	case i64:
		dst = i64(i)
	case u8:
		dst = u8(i)
	case u16:
		dst = u16(i)
	case u32:
		dst = u32(i)
	case u64:
		dst = u64(i)
	case int:
		dst = int(i)
	case uint:
		dst = uint(i)
	case:
		assert(false, "tried to assign int on something that cant be assigned int")
	}
}

@(private)
assign_float :: proc(dest: any, f: f64) {
	v := reflect.any_core(dest)
	switch &dst in v {
	case f32:
		dst = f32(f)
	case f64:
		dst = f64(f)
	case:
		assert(false, "tried to assign float on something that cant be assigned float")
	}
}

@(private)
unmarshal_value :: proc(dest: any, value: toml.Type, path: string) -> bool {
	ti := reflect.type_info_base(type_info_of(dest.id))

	#partial switch t in ti.variant {
	case reflect.Type_Info_Boolean:
		v, ok := value.(bool)
		if !ok {
			log.errorf("config: %q expected bool, got %T", path, value)
			return false
		}
		(^bool)(dest.data)^ = v

	case reflect.Type_Info_Integer:
		v, ok := value.(i64)
		if !ok {
			log.errorf("config: %q expected integer, got %T", path, value)
			return false
		}
		if !t.signed && v < 0 {
			log.errorf("config: %q expected positive integer, got %v", path, v)
			return false
		}
		assign_int(dest, v)

	case reflect.Type_Info_Float:
		#partial switch v in value {
		case f64:
			assign_float(dest, v)
		case i64:
			assign_float(dest, f64(v))
		case:
			log.errorf("config: %q expected float, got %T", path, value)
			return false
		}

	case reflect.Type_Info_String:
		v, ok := value.(string)
		if !ok {
			log.errorf("config: %q expected string, got %T", path, value)
			return false
		}
		(^string)(dest.data)^ = v

	case reflect.Type_Info_Struct:
		v, ok := value.(^toml.Table)
		if !ok {
			log.errorf("config: %q expected a table, got %T", path, value)
			return false
		}
		return unmarshal_table(dest, v, path)

	case reflect.Type_Info_Pointer:
		if t.elem.id == typeid_of(toml.Table) {
			v, ok := value.(^toml.Table)
			if !ok {
				log.errorf("config: %q expected table, got %T", path, value)
				return false
			}
			(^(^toml.Table))(dest.data)^ = v
		} else if t.elem.id == typeid_of(toml.List) {
			v, ok := value.(^toml.List)
			if !ok {
				log.errorf("config: %q expected list, got %T", path, value)
				return false
			}
			(^(^toml.List))(dest.data)^ = v
		} else {
			assert(
				false,
				"config: unsupported pointer type in schema — only ^toml.Table and ^toml.List are allowed, use value types for nested structs",
			)
		}
	case reflect.Type_Info_Enum:
		v, ok := value.(string)
		if !ok {
			log.errorf("config: %q expected string for enum, got %T", path, value)
			return false
		}
		for name, i in t.names {
			if name == v {
				assign_int(dest, i64(t.values[i]))
				return true
			}
		}
		log.errorf("config: %q invalid enum value %q expected one of %v", path, v, t.names)
		return false

	case:
		assert(false, "unsupported field type in schema")
	}

	return true
}

@(private)
unmarshal_table :: proc(dest: any, table: ^toml.Table, path: string) -> bool {
	ti := reflect.type_info_base(type_info_of(dest.id))
	t, ok := ti.variant.(reflect.Type_Info_Struct)
	assert(ok, "unmarshal destination must be a struct")
	assert(.raw_union not_in t.flags, "raw_union structs are not supported")

	fields := reflect.struct_fields_zipped(ti.id)

	for field in fields {
		tag_val := reflect.struct_tag_get(field.tag, "toml")
		key, required := parse_tag(tag_val)
		if key == "" {key = field.name}

		field_path := fmt.tprintf("%s.%s", path, key)

		value, exists := table[key]
		if !exists {
			if required {
				log.errorf("config: required field %q is missing", field_path)
				return false
			}
			continue
		}

		field_ptr := rawptr(uintptr(dest.data) + field.offset)
		field_any := any {
			data = field_ptr,
			id   = field.type.id,
		}
		if !unmarshal_value(field_any, value, field_path) {
			return false
		}
	}

	return true
}

@(private)
table_get_raw :: proc(
	$T: typeid,
	section: ^toml.Table,
	optional: bool,
	path: ..string,
) -> (
	val: T,
	ok: bool,
) where intrinsics.type_is_variant_of(toml.Type, T) {
	assert(len(path) > 0, "You must specify at least one path str in table_get()!")


	full_path := strings.join(path, ".", context.temp_allocator)

	if section == nil {
		if !optional {log.errorf("config: %q table is nil", full_path)}
		return {}, false
	}

	section := section
	for dir in path[:len(path) - 1] {
		next, next_ok := section[dir].(^toml.Table)
		if !next_ok {
			if !optional {log.errorf("config: field %q does not exist", full_path)}
			return {}, false
		}
		section = next
	}

	last := path[len(path) - 1]
	raw, exists := section[last]
	if !exists {
		if !optional {log.errorf("config: field %q does not exist", full_path)}
		return {}, false
	}

	v, val_ok := raw.(T)
	if !val_ok {
		if !optional {log.errorf("config: %q expected %v, got %T", full_path, typeid_of(T), raw)}
		return {}, false
	}
	val = v

	return val, true
}

table_get :: proc($T: typeid, section: ^toml.Table, path: ..string) -> (val: T, ok: bool) {
	return table_get_raw(T, section, false, ..path)
}
table_get_opt :: proc($T: typeid, section: ^toml.Table, path: ..string) -> (val: T, ok: bool) {
	return table_get_raw(T, section, true, ..path)
}


unmarshal :: proc(
	table: ^toml.Table,
	ptr: ^$T,
	temp_allocator := context.temp_allocator,
) -> bool where intrinsics.type_is_struct(T) {
	context.temp_allocator = temp_allocator
	assert(table != nil, "config: unmarshal called with nil table")
	return unmarshal_table(any{data = ptr, id = typeid_of(T)}, table, fmt.tprintf("%T", ptr^))
}
