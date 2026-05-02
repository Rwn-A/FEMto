package fem

import "core:math"
import "core:mem"
import "core:slice"

Time_Scheme :: enum {
	Steady,
	BE,
	BDF2,
	Newmark,
}

Timestepper :: struct {
	start, end, dt, current_time: f64,
	current_step, total_steps:    int,
	schemes:                      [MAX_VARS]Time_Scheme,
	allocator:                    mem.Allocator,
}

Timestep_State :: struct {
	previous_soln:          Vector,
	previous_previous_soln: Vector,
	previous_ddt:           Vector,
	previous_d2dt2:         Vector,
}

Timestep :: struct {
	step:      int,
	dt:        f64,
	time:      f64,
	du_du_dot: [MAX_VARS]f64,
}

timestepper_create :: proc(
	start, end, dt: f64,
	ics: Vector,
	allocator := context.allocator,
) -> (
	Timestepper,
	Timestep_State,
) {
	ts := Timestepper {
		start        = start,
		end          = end,
		dt           = dt,
		current_time = start,
		total_steps  = int(math.ceil((end - start) / dt)),
		allocator    = allocator,
	}
	state := Timestep_State {
		previous_soln = slice.clone(ics, allocator),
	}
	return ts, state
}

// convience layer, sets all schemes to steady, only does on step, all time derivative quantities will be 0.
// State is still allocated with the room for the previous solution.
timestepper_create_steady :: proc(
	ics: Vector,
	sys: System,
	allocator := context.allocator,
) -> (
	Timestepper,
	Timestep_State,
) {
	ts := Timestepper {
		dt          = 1.0,
		total_steps = 1,
		allocator   = allocator,
	}
	state := Timestep_State {
		previous_soln = slice.clone(ics, allocator),
	}

	slice.fill(ts.schemes[:], Time_Scheme.Steady)

	return ts, state
}

timestep_state_destroy :: proc(ts: ^Timestepper, state: Timestep_State) {
	context.allocator = ts.allocator
	delete(state.previous_soln)
	if state.previous_previous_soln != nil {delete(state.previous_previous_soln)}
	if state.previous_ddt != nil {delete(state.previous_ddt)}
	if state.previous_d2dt2 != nil {delete(state.previous_d2dt2)}
}

timestepper_set_scheme :: proc(ts: ^Timestepper, state: ^Timestep_State, handle: Var_Handle, scheme: Time_Scheme) {
	ts.schemes[handle] = scheme
	n := len(state.previous_soln)
	switch scheme {
	case .BDF2:
		if state.previous_previous_soln == nil {
			state.previous_previous_soln = make(Vector, n, ts.allocator)
		}
	case .Newmark:
		if state.previous_ddt == nil {
			state.previous_ddt = make(Vector, n, ts.allocator)
		}
		if state.previous_d2dt2 == nil {
			state.previous_d2dt2 = make(Vector, n, ts.allocator)
		}
	case .BE, .Steady:
	// previous_soln already allocated in create
	}
}

timestepper_du_du_dot :: proc(ts: ^Timestepper, step: Timestep, handle: Var_Handle) -> f64 {
	switch ts.schemes[handle] {
	case .Steady:
		return 0.0
	case .BE:
		return 1.0 / step.dt
	case .BDF2:
		if ts.current_step <= 1 {
			return 1.0 / step.dt
		}
		return 3.0 / (2.0 * step.dt)
	case .Newmark:
		beta :: 0.25
		return 1.0 / (beta * step.dt * step.dt)
	}
	unreachable()
}

timestepper_step :: proc(ts: ^Timestepper, state: Timestep_State, u: Vector, sys: System) -> (Timestep, bool) {
	if ts.current_step >= ts.total_steps {return {}, false}
	if ts.current_step > 0 {
		u_dot, u_ddot := timestepper_derivatives(ts, state, u, sys, ts.allocator)
		defer {delete(u_dot, ts.allocator); delete(u_ddot, ts.allocator)}

		if state.previous_previous_soln != nil {
			copy(state.previous_previous_soln, state.previous_soln)
		}

		copy(state.previous_soln, u)

		if state.previous_ddt != nil {
			copy(state.previous_ddt, u_dot)
		}
		if state.previous_d2dt2 != nil {
			copy(state.previous_d2dt2, u_ddot)
		}
	}
	step := Timestep {
		time = ts.current_time,
		step = ts.current_step,
		dt   = ts.dt,
	}
	for var in 0 ..< sys.num_vars {
		step.du_du_dot[var] = timestepper_du_du_dot(ts, step, Var_Handle(var))
	}
	defer {ts.current_step += 1; ts.current_time += ts.dt}
	return step, true
}

timestepper_derivatives :: proc(
	ts: ^Timestepper,
	state: Timestep_State,
	u: Vector,
	sys: System,
	allocator := context.allocator,
) -> (
	u_dot, u_ddot: Vector,
) {
	u_dot = make(Vector, len(u), allocator)
	u_ddot = make(Vector, len(u), allocator)
	for var in 0 ..< sys.num_vars {
		u_cur_range := system_var_slice(sys, var, u)
		u_dot_range := system_var_slice(sys, var, u_dot)
		u_ddot_range := system_var_slice(sys, var, u_ddot)
		u_prev_range := system_var_slice(sys, var, state.previous_soln)
		switch ts.schemes[var] {
		case .Steady:
		case .BE:
			copy(u_dot_range, u_cur_range)
			scal(u_dot_range, 1.0 / ts.dt)
			axpy(u_prev_range, u_dot_range, -1.0 / ts.dt)
		case .BDF2:
			if ts.current_step <= 1 || state.previous_previous_soln == nil {
				copy(u_dot_range, u_cur_range)
				scal(u_dot_range, 1.0 / ts.dt)
				axpy(u_prev_range, u_dot_range, -1.0 / ts.dt)
			} else {
				u_prev_prev_range := system_var_slice(sys, var, state.previous_previous_soln)
				copy(u_dot_range, u_cur_range)
				scal(u_dot_range, 3.0 / (2.0 * ts.dt))
				axpy(u_prev_range, u_dot_range, -4.0 / (2.0 * ts.dt))
				axpy(u_prev_prev_range, u_dot_range, 1.0 / (2.0 * ts.dt))
			}
		case .Newmark:
			beta :: 0.25
			gamma :: 0.5
			u_dot_prev_range := system_var_slice(sys, var, state.previous_ddt)
			u_ddot_prev_range := system_var_slice(sys, var, state.previous_d2dt2)
			for i in 0 ..< len(u_dot_range) {
				u_ddot_range[i] =
					(u_cur_range[i] - u_prev_range[i]) / (beta * ts.dt * ts.dt) -
					u_dot_prev_range[i] / (beta * ts.dt) -
					(1.0 - 2.0 * beta) / (2.0 * beta) * u_ddot_prev_range[i]
				u_dot_range[i] =
					u_dot_prev_range[i] +
					ts.dt * (1.0 - gamma) * u_ddot_prev_range[i] +
					ts.dt * gamma * u_ddot_range[i]
			}
		}
	}
	return u_dot, u_ddot
}
