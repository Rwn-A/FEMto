/*
 Wraps the vendor/toml library to support better error messages and richer unmarshalling.

 USAGE:
 ```
      cfg := struct {
          host_xyz: string `toml:"host" validate:"required"`,
          port: int        `validate:"gte=0,lte=65535"`,
          tags: []string,
      } {
        // default values can go here if needed.
        port = 8080
      }

      ok := unmarshal(table, &cfg)
```

NOTES:
- Strings are cloned by the unmarshaller, while ^Table and ^Lists are not.
- The range validation (lt, lte, gt, gte) can only be applied to integers, floats and arrays of them.
- Invalid ranges are not dected.  (ex. lte=10,gte=12)

*/
package cfg

import "base:intrinsics"
import "core:fmt"
import "core:log"
import "core:mem"
import "core:reflect"
import "core:strconv"
import "core:strings"

import "../../vendor/toml"

// In the case other code interacts with these toml layer doesnt leak.
Table :: toml.Table
List :: toml.List
Type :: toml.Type
parse :: toml.parse
format_error :: toml.format_error

// Typechecks and retrieves the value specified by `path`.
// Logs an error and return false if not found.
// Note: uses context.temp_allocator in the case of an error.
table_get :: proc($T: typeid, section: ^Table, path: ..string) -> (T, bool) {
	context.allocator = context.temp_allocator
	return table_get_raw(T, section, false, ..path)
}

// Typechecks and retrieves the value specified by `path`, returns false is the value doesnt exist
// this function still logs an error if the value exists but does not match the expected type.
table_get_opt :: proc($T: typeid, section: ^Table, path: ..string) -> (T, bool) {
	context.allocator = context.temp_allocator
	return table_get_raw(T, section, true, ..path)
}

// Tries to load fields from configuration in `table` into the struct pointed to by `ptr`.
// All memory is allocated with the allocator, there is no delete proc, it is expected to be freed at once.
// If a `Table` or `List` is left as is through the schema it is not copied.
// The schema struct can optionally use struct tags for further verification (see usage above).
// Error messages are automatically logged, the returned boolean indicates if this happened.
// If the schema struct has fields that cannot be unmarshalled an assert will be triggered.
unmarshal :: proc(table: ^Table, ptr: ^$T, allocator := context.allocator) -> bool where intrinsics.type_is_struct(T) {
	assert(table != nil, fmt.aprintf("unmarshal called with nil table %T", ptr^))
	context.allocator = allocator
	return unmarshal_table(any{data = ptr, id = typeid_of(T)}, table, fmt.aprintf("%T", ptr^))
}


@(private = "file")
Validation_Options :: bit_set[Validation_Option]

@(private = "file")
NUMERIC_VALIDATORS :: Validation_Options{.GT, .GTE, .LT, .LTE}


@(private = "file")
Validation_Option :: enum {
	Required,
	GT,
	LT,
	GTE,
	LTE,
}

@(private = "file")
Validator :: struct {
	opt:              Validation_Options,
	gt, gte, lt, lte: f64,
}

@(private = "file")
parse_toml_tag :: proc(tag: reflect.Struct_Tag) -> string {
	v := reflect.struct_tag_get(tag, "toml")
	return v
}

@(private = "file")
parse_validate_tag :: proc(tag: reflect.Struct_Tag) -> (v: Validator) {
	raw := reflect.struct_tag_get(tag, "validate")
	if raw == "" {return {}}

	parts := strings.split(raw, ",")
	for part in parts {
		part := strings.trim_space(part)
		if part == "required" {
			v.opt += {.Required}
			continue
		}
		eq := strings.index_byte(part, '=')
		assert(eq != -1, fmt.aprintf("config: unrecognised validate token %q", part))

		key := part[:eq]
		val_str := part[eq + 1:]
		val, parse_ok := strconv.parse_f64(val_str)
		assert(parse_ok, fmt.aprintf("config: validate token %q has non-numeric value %q", key, val_str))

		switch key {
		case "gt":
			v.opt += {.GT}; v.gt = val
		case "gte":
			v.opt += {.GTE}; v.gte = val
		case "lt":
			v.opt += {.LT}; v.lt = val
		case "lte":
			v.opt += {.LTE}; v.lte = val
		case:
			assert(false, fmt.aprintf("config: unrecognised validate key %q", key))
		}
	}

	return v
}

@(private = "file")
validate_number :: proc(vd: Validator, num: f64, path: string) -> bool {
	if .GT in vd.opt && !(num > vd.gt) {log.errorf("config: %q value %v must be > %v", path, num, vd.gt); return false}
	if .GTE in vd.opt &&
	   !(num >= vd.gte) {log.errorf("config: %q value %v must be >= %v", path, num, vd.gte); return false}
	if .LT in vd.opt && !(num < vd.lt) {log.errorf("config: %q value %v must be < %v", path, num, vd.lt); return false}
	if .LTE in vd.opt &&
	   !(num <= vd.lte) {log.errorf("config: %q value %v must be <= %v", path, num, vd.lte); return false}
	return true
}


// Below unmarshalling based on original unmarshalling in the toml library as well as core:encoding/json.

@(private = "file")
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
		assert(false, "tried to assign int to an unsupported type")
	}
}

@(private = "file")
assign_float :: proc(dest: any, f: f64) {
	v := reflect.any_core(dest)
	switch &dst in v {
	case f32:
		dst = f32(f)
	case f64:
		dst = f64(f)
	case:
		assert(false, "tried to assign float to an unsupported type")
	}
}

@(private = "file")
numeric_as_f64 :: proc(value: toml.Type) -> (f64, bool) {
	#partial switch v in value {
	case f64:
		return v, true
	case i64:
		return f64(v), true
	case:
		return 0, false
	}
}

@(private)
unmarshal_value :: proc(dest: any, value: toml.Type, path: string, vd: Validator) -> bool {
	ti := reflect.type_info_base(type_info_of(dest.id))

	#partial switch t in ti.variant {
	case reflect.Type_Info_Boolean:
		v, ok := value.(bool)
		if !ok {log.errorf("config: %q expected bool, got %T", path, value); return false}
		(^bool)(dest.data)^ = v

	case reflect.Type_Info_Integer:
		v, ok := value.(i64)
		if !ok {log.errorf("config: %q expected integer, got %T", path, value); return false}
		if !t.signed && v < 0 {
			log.errorf("config: %q expected positive integer, got %v", path, v)
			return false
		}
		if !validate_number(vd, f64(v), path) {return false}
		assign_int(dest, v)

	case reflect.Type_Info_Float:
		num, ok := numeric_as_f64(value)
		if !ok {log.errorf("config: %q expected float, got %T", path, value); return false}
		if !validate_number(vd, num, path) {return false}
		assign_float(dest, num)

	case reflect.Type_Info_String:
		v, ok := value.(string)
		if !ok {log.errorf("config: %q expected string, got %T", path, value); return false}
		(^string)(dest.data)^ = strings.clone(v)

	case reflect.Type_Info_Struct:
		v, ok := value.(^Table)
		if !ok {log.errorf("config: %q expected a table, got %T", path, value); return false}
		return unmarshal_table(dest, v, path)

	case reflect.Type_Info_Enum:
		v, ok := value.(string)
		if !ok {log.errorf("config: %q expected string for enum, got %T", path, value); return false}
		for name, i in t.names {
			if name == v {
				assign_int(dest, i64(t.values[i]))
				return true
			}
		}
		log.errorf("config: %q invalid enum value %q, expected one of %v", path, v, t.names)
		return false

	case reflect.Type_Info_Array:
		list, ok := value.(^List)
		if !ok {log.errorf("config: %q expected a list, got %T", path, value); return false}
		if len(list) != t.count {
			log.errorf("config: %q array length mismatch: expected %d elements, got %d", path, t.count, len(list))
			return false
		}
		return assign_list(dest.data, t.elem, list, path, vd)

	case reflect.Type_Info_Slice:
		list, ok := value.(^List)
		if !ok {log.errorf("config: %q expected a list, got %T", path, value); return false}
		raw := (^mem.Raw_Slice)(dest.data)
		data, data_ok := mem.alloc_bytes(t.elem.size * len(list), t.elem.align)
		if data_ok != .None {
			log.errorf("config: %q failed to allocate slice (%v)", path, data_ok)
			return false
		}
		raw.data = raw_data(data)
		raw.len = len(list)
		return assign_list(raw.data, t.elem, list, path, vd)

	// Allows List and Table to be left intact for more dynamic parsing using `table_get` or similar.
	case reflect.Type_Info_Pointer:
		if t.elem.id == typeid_of(Table) {
			v, ok := value.(^Table)
			if !ok {log.errorf("config: %q expected table, got %T", path, value); return false}
			(^(^Table))(dest.data)^ = v
		} else if t.elem.id == typeid_of(List) {
			v, ok := value.(^List)
			if !ok {log.errorf("config: %q expected list, got %T", path, value); return false}
			(^(^List))(dest.data)^ = v
		} else {
			assert(
				false,
				"config: unsupported pointer type — only ^Table and ^List are allowed; use value types for nested structs",
			)

		}
	case:
		assert(false, "config: unsupported field type in schema")
	}

	return true
}

@(private = "file")
assign_list :: proc(base: rawptr, elem: ^reflect.Type_Info, list: ^List, path: string, vd: Validator) -> bool {
	for i in 0 ..< len(list) {
		elem_path := fmt.aprintf("%s[%d]", path, i)
		elem_ptr := rawptr(uintptr(base) + uintptr(i) * uintptr(elem.size))
		elem_any := any {
			data = elem_ptr,
			id   = elem.id,
		}
		if !unmarshal_value(elem_any, list[i], elem_path, vd) {return false}
	}
	return true
}

@(private = "file")
unmarshal_table :: proc(dest: any, table: ^Table, path: string) -> bool {
	ti := reflect.type_info_base(type_info_of(dest.id))
	t, ok := ti.variant.(reflect.Type_Info_Struct)
	assert(ok, "unmarshal destination must be a struct")
	assert(.raw_union not_in t.flags, "raw_union structs are not supported")

	for field in reflect.struct_fields_zipped(ti.id) {
		key := parse_toml_tag(field.tag)
		if key == "" {key = field.name}

		field_path := fmt.aprintf("%s.%s", path, key)

		vd := parse_validate_tag(field.tag)

		field_ti := reflect.type_info_base(field.type)
		_, is_int := field_ti.variant.(reflect.Type_Info_Integer)
		_, is_float := field_ti.variant.(reflect.Type_Info_Float)
		_, is_arr := field_ti.variant.(reflect.Type_Info_Array)
		_, is_slice := field_ti.variant.(reflect.Type_Info_Slice)
		_, is_enum := field_ti.variant.(reflect.Type_Info_Enum)
		numeric_validators_present := vd.opt & NUMERIC_VALIDATORS != {}
		if numeric_validators_present {
			assert(
				is_int || is_float || is_arr || is_slice,
				fmt.aprintf(
					"config: numeric validators on %q are only valid for int, float, array, or slice fields",
					field_path,
				),
			)
			if is_enum {
				assert(false, fmt.aprintf("config: numeric validators not allowed on enum field %q", field_path))
			}
		}

		value, exists := table[key]
		if !exists {
			if .Required in vd.opt {
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
		if !unmarshal_value(field_any, value, field_path, vd) {
			return false
		}
	}

	return true
}


// Based on the `get` function from the the toml library.
@(private = "file")
table_get_raw :: proc(
	$T: typeid,
	section: ^Table,
	optional: bool,
	path: ..string,
) -> (
	val: T,
	ok: bool,
) where intrinsics.type_is_variant_of(toml.Type, T) {
	assert(len(path) > 0, "table_get: must provide at least one path segment")

	if section == nil {
		if !optional {
			full_path := strings.join(path, ".")
			log.errorf("config: %q table is nil", full_path)
		}
		return {}, false
	}

	section := section
	for dir in path[:len(path) - 1] {
		next, next_ok := section[dir].(^Table)
		if !next_ok {
			if !optional {
				full_path := strings.join(path, ".")
				log.errorf("config: field %q does not exist", full_path)
			}
			return {}, false
		}
		section = next
	}

	last := path[len(path) - 1]
	raw, exists := section[last]
	if !exists {
		if !optional {
			full_path := strings.join(path, ".")
			log.errorf("config: field %q does not exist", full_path)
		}
		return {}, false
	}

	v, val_ok := raw.(T)
	if !val_ok {
		full_path := strings.join(path, ".")
		log.errorf("config: %q expected %v, got %T", full_path, typeid_of(T), raw)
		return {}, false
	}
	return v, true
}
