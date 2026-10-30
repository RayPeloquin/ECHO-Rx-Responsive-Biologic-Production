"""Reproduce the reduced ECHO-Rx shutdown illustration with Python 3.9+.

Added during the 29 October 2026 repository repair from the published equations.
All masses and times are normalized, arbitrary examples. This is not a biological
simulation, dose calculation, clearance guarantee, or model of residual secretion.
Run with no arguments to print the checked JSON record. No files are overwritten.
"""
from dataclasses import asdict, dataclass
import json
import math


@dataclass(frozen=True)
class Parameters:
    depot_initial: float = 6.0
    target_initial: float = 1.0
    a: float = 0.2
    b: float = 0.15
    c: float = 0.5

    def __post_init__(self):
        values = asdict(self).values()
        if any(type(v) not in (int, float) or not math.isfinite(v) for v in values):
            raise ValueError("Parameters must be finite real numbers")
        if (self.depot_initial < 0 or self.target_initial < 0
                or self.a <= 0 or self.c <= 0 or not 0 <= self.b <= self.a):
            raise ValueError("Require nonnegative masses, a,c > 0 and 0 <= b <= a")


def analytical(time, p):
    """Return depot and target mass after an ideal instantaneous secretion stop."""
    if type(time) not in (int, float) or not math.isfinite(time) or time < 0:
        raise ValueError("Time must be finite and nonnegative")
    depot = p.depot_initial * math.exp(-p.a * time)
    difference = abs(p.c - p.a)
    if difference == 0:
        transfer = time * math.exp(-p.a * time)
    else:
        # Stable divided difference, including almost-equal rates. The exponent
        # is never positive, avoiding cancellation and exponential overflow.
        transfer = (math.exp(-min(p.a, p.c) * time)
                    * (-math.expm1(-difference * time)) / difference)
    target = (p.target_initial * math.exp(-p.c * time)
              + p.b * p.depot_initial * transfer)
    return depot, target


def remaining_target_integral(depot, target, p):
    """Integral of target mass from this state to infinity, under this model."""
    return (target + p.b * depot / p.a) / p.c


def rk4(p, horizon=80.0, step=0.01):
    """Independent ODE integration; states are depot, target, loss and integral."""
    if not (math.isfinite(horizon) and math.isfinite(step) and horizon > 0 and step > 0):
        raise ValueError("Positive finite horizon and step required")
    count = math.ceil(horizon / step)
    step = horizon / count
    state = (p.depot_initial, p.target_initial, 0.0, 0.0)
    records = [(0.0, state)]

    def derivative(y):
        depot, target, _, _ = y
        return (-p.a * depot, p.b * depot - p.c * target,
                (p.a - p.b) * depot + p.c * target, target)

    def advance(y, change, factor):
        return tuple(v + factor * dv for v, dv in zip(y, change))

    for i in range(count):
        k1 = derivative(state)
        k2 = derivative(advance(state, k1, step / 2))
        k3 = derivative(advance(state, k2, step / 2))
        k4 = derivative(advance(state, k3, step))
        state = tuple(v + step * (a + 2*b + 2*c + d) / 6
                      for v, a, b, c, d in zip(state, k1, k2, k3, k4))
        records.append(((i + 1) * step, state))
    return records


def checked_record():
    checks = []

    def require(condition, label):
        # Explicit conditions remain active under python -O.
        if not condition:
            raise RuntimeError("Check failed: " + label)
        checks.append(label)

    p = Parameters()
    records = rk4(p)
    max_error = max(max(abs(y[i] - analytical(t, p)[i]) for i in (0, 1))
                    for t, y in records)
    require(max_error < 1e-9, "Independent RK4 matches the analytical trajectory")
    total = p.depot_initial + p.target_initial
    balance_error = max(abs(sum(y[:3]) - total) for _, y in records)
    require(balance_error < 1e-10, "Depot plus target plus eliminated mass is conserved")
    peak_time = math.log(5.0 / 3.0) / 0.3
    peak_depot, peak_target = analytical(peak_time, p)
    require(abs(p.b * peak_depot - p.c * peak_target) < 1e-12
            and abs(peak_target - 1.2804815896) < 1e-10,
            "The documented post-stop peak and its zero slope agree")
    integral = remaining_target_integral(p.depot_initial, p.target_initial, p)
    _, final = records[-1]
    recovered_integral = final[3] + remaining_target_integral(final[0], final[1], p)
    require(abs(integral - 11.0) < 1e-12 and abs(recovered_integral - integral) < 1e-9,
            "Integrated exposure including the remaining tail equals 11")
    equal = Parameters(a=0.2, b=0.15, c=0.2)
    equal_error = max(abs(y[1] - analytical(t, equal)[1]) for t, y in rk4(equal))
    near = Parameters(a=0.2, b=0.15, c=0.2000000000001)
    require(equal_error < 1e-9 and all(
        abs(analytical(t, equal)[1] - analytical(t, near)[1]) < 1e-10
        for t in (0, 1, 5, 20)), "Equal and nearly equal rates have a stable limit")
    steady = Parameters(target_initial=p.b * p.depot_initial / p.c)
    trajectory = rk4(steady)
    require(all(after[1][1] <= before[1][1] + 1e-12
                for before, after in zip(trajectory, trajectory[1:])),
            "An exact pre-OFF steady state decays without a later target peak")
    zero = Parameters(depot_initial=0.0, b=0.0)
    require(all(abs(analytical(t, zero)[1] - math.exp(-p.c*t)) < 1e-14
                for t in (0, 1, 10)), "No-depot transfer reduces to simple target decay")
    for invalid in ({"a": 0}, {"c": -1}, {"b": 0.3},
                    {"depot_initial": -1}, {"target_initial": float("nan")}):
        try:
            Parameters(**invalid)
        except ValueError:
            continue
        raise RuntimeError("Invalid parameters were accepted")
    checks.append("Nonfinite or physically inconsistent reduced-model parameters are rejected")
    return {
        "schema": "echo-rx-reduced-shutdown/1",
        "status": "Illustrative mathematics only; no biological or clinical measurements",
        "units": "Normalized mass and arbitrary time",
        "parameters": asdict(p),
        "example": {
            "target_expression": "3 exp(-0.2t) - 2 exp(-0.5t)",
            "initial_target_slope": p.b*p.depot_initial - p.c*p.target_initial,
            "peak_time": peak_time, "peak_target": peak_target,
            "target_integral_to_infinity": integral,
        },
        "numerical_check": {
            "method": "RK4; depot, target, eliminated mass and target integral",
            "horizon": 80.0, "step": 0.01,
            "maximum_absolute_trajectory_error": max_error,
            "maximum_mass_balance_error": balance_error,
            "equal_rate_maximum_absolute_target_error": equal_error,
        },
        "checks": checks,
        "limitations": [
            "Ideal instantaneous stop of all further secretion; no leak or death release",
            "Fixed linear coefficients and homogeneous compartments",
            "Only the two-compartment illustration, not the full biological model",
            "No dose, biological clearance time or patient-use interpretation",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(checked_record(), indent=2, allow_nan=False))
