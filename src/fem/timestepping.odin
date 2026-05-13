package fem

import "core:math"
import "core:mem"
import "core:slice"

// TODO: more control over newmark, velocity ic's,

NEWMARK_GAMMA :: 0.5
NEWMARK_BETA :: 0.25

Timestepper :: struct {
	start, end, dt, current_time: f64,
	current_step, total_steps:    int,
	schemes:                      [MAX_VARS]Time_Scheme,
	history: Timestepper_History,
}

Timestepper_History :: struct {
	previous_soln:          Vector,
	previous_previous_soln: Vector,
	previous_ddt:           Vector,
	previous_d2dt2:         Vector,
}

Timestep_State :: struct {
	step: Timestep,
	ddt, d2dt2: Vector,
}

Time_Scheme :: enum {
	Steady,
	BE,
	BDF2,
	Newmark,
}

Timestep :: struct {
	step:       int,
	dt:         f64,
	time:       f64,
	du_du_dot:  [MAX_VARS]f64,
	du_du_ddot: [MAX_VARS]f64,
}

timestepper_create :: proc(
	start, end, dt: f64,
) -> (
	Timestepper,
) {
	ts := Timestepper {
		start        = start,
		end          = end,
		dt           = dt,
		current_time = start + dt,
		current_step = 1,
		total_steps  = int(math.ceil((end - start) / dt)),
	}
	return ts
}

// convience layer, sets all schemes to steady, only does on step, all time derivative quantities will be 0.
timestepper_create_steady :: proc(
	allocator := context.allocator,
) -> (
	Timestepper,
) {
	ts := Timestepper {
		dt           = 1.0,
		total_steps  = 1,
		current_step = 1,
	}

	slice.fill(ts.schemes[:], Time_Scheme.Steady)

	return ts
}


timestepper_set_scheme :: proc(ts: ^Timestepper, handle: Var_Handle, scheme: Time_Scheme) {
	ts.schemes[handle] = scheme
}

// used for IC's, does not contain all the fields of a timestep, just `step` and `time`.
timestepper_initial_step :: proc(ts: Timestepper) -> Timestep {
	return {step = 0, time = ts.start}
}

// must be called after your schemes have beens set.
timestepper_init_history :: proc(ts: ^Timestepper, sys: System, ics: Vector, allocator := context.allocator) {
	context.allocator = allocator

	ts := ts
	schemes := slice.enum_slice_to_bitset(ts.schemes[0:sys.num_vars], bit_set[Time_Scheme])

	ts.history.previous_soln = slice.clone(ics)

	if .BDF2 in schemes {
		ts.history.previous_previous_soln = system_vector(sys)
	}

	if .Newmark in schemes {
		ts.history.previous_ddt = system_vector(sys)
		ts.history.previous_d2dt2 = system_vector(sys)
	}

}

timestep_destroy_history :: proc(ts: ^Timestepper, allocator := context.allocator) {
	context.allocator = allocator
	delete(ts.history.previous_soln)
	if ts.history.previous_previous_soln != nil {delete(ts.history.previous_previous_soln)}
	if ts.history.previous_ddt != nil {delete(ts.history.previous_ddt)}
	if ts.history.previous_d2dt2 != nil {delete(ts.history.previous_d2dt2)}
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
		return NEWMARK_GAMMA / (NEWMARK_BETA * step.dt)
	}
	unreachable()
}


timestepper_du_du_ddot :: proc(ts: ^Timestepper, step: Timestep, handle: Var_Handle) -> f64 {
	switch ts.schemes[handle] {
	case .Steady, .BE, .BDF2:
		return 0
	case .Newmark:
		return 1.0 / (NEWMARK_BETA * step.dt * step.dt)
	}
	unreachable()
}

timestepper_step :: proc(ts: ^Timestepper, sys: System, u: Vector, step_state: ^Timestep_State) -> (Timestep, bool) {
	if ts.current_step > ts.total_steps {return {}, false}
	if ts.current_step > 1 {
		timestepper_derivatives(ts, sys, u, step_state.ddt, step_state.d2dt2)

		if ts.history.previous_previous_soln != nil {
			copy(ts.history.previous_previous_soln, ts.history.previous_soln)
		}

		copy(ts.history.previous_soln, u)

		if ts.history.previous_ddt != nil {
			copy(ts.history.previous_ddt, step_state.ddt)
		}
		if ts.history.previous_d2dt2 != nil {
			copy(ts.history.previous_d2dt2, step_state.d2dt2)
		}
	}
	step_state.step = Timestep {
		time = ts.current_time,
		step = ts.current_step,
		dt   = ts.dt,
	}
	for var in 0 ..< sys.num_vars {
		step_state.step.du_du_dot[var] = timestepper_du_du_dot(ts, step_state.step, Var_Handle(var))
		step_state.step.du_du_ddot[var] = timestepper_du_du_ddot(ts, step_state.step, Var_Handle(var))
	}

	defer {ts.current_step += 1; ts.current_time += ts.dt}

	return step_state.step, true
}

// Required to be called if your solution changes within a timestep (coupling / nonlinear iterations)
timestepper_update_derivatives :: proc(ts: ^Timestepper, step_state: ^Timestep_State, sys: System, u: Vector) {
	timestepper_derivatives(ts, sys, u, step_state.ddt, step_state.d2dt2)
}

timestepper_derivatives :: proc(
	ts: ^Timestepper,
	sys: System,
	u: Vector,
	u_dot, u_ddot: Vector,
) {
	for var in 0 ..< sys.num_vars {
		u_cur_range := system_var_slice(sys, var, u)
		u_dot_range := system_var_slice(sys, var, u_dot)
		u_ddot_range := system_var_slice(sys, var, u_ddot)
		u_prev_range := system_var_slice(sys, var, ts.history.previous_soln)
		switch ts.schemes[var] {
		case .Steady:
		case .BE:
			copy(u_dot_range, u_cur_range)
			scal(u_dot_range, 1.0 / ts.dt)
			axpy(u_prev_range, u_dot_range, -1.0 / ts.dt)
		case .BDF2:
			if ts.current_step <= 1 || ts.history.previous_previous_soln == nil {
				copy(u_dot_range, u_cur_range)
				scal(u_dot_range, 1.0 / ts.dt)
				axpy(u_prev_range, u_dot_range, -1.0 / ts.dt)
			} else {
				u_prev_prev_range := system_var_slice(sys, var, ts.history.previous_previous_soln)
				copy(u_dot_range, u_cur_range)
				scal(u_dot_range, 3.0 / (2.0 * ts.dt))
				axpy(u_prev_range, u_dot_range, -4.0 / (2.0 * ts.dt))
				axpy(u_prev_prev_range, u_dot_range, 1.0 / (2.0 * ts.dt))
			}
		case .Newmark:
			beta :: 0.25
			gamma :: 0.5
			u_dot_prev_range := system_var_slice(sys, var, ts.history.previous_ddt)
			u_ddot_prev_range := system_var_slice(sys, var, ts.history.previous_d2dt2)
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
}
