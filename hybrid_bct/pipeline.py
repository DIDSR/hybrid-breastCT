from .systems.doheny import DohenySystem


def run_hybrid_simulation(
    cfg: dict,
    scan_id: int,
    cluster_diameter_mm: float,
    num_calcs: int,
    calc_diameter_mm: float,
) -> None:
    system_name = cfg.get("system", "").lower()

    if system_name == "doheny":
        system = DohenySystem(cfg)
    else:
        raise ValueError(f"Unsupported system: {system_name}")

    system.validate()

    print("Configuration validated.")
    print(f"System summary: {system.summary()}")
    print(
        f"Run request: scan_id={scan_id}, "
        f"cluster_diameter_mm={cluster_diameter_mm}, "
        f"num_calcs={num_calcs}, "
        f"calc_diameter_mm={calc_diameter_mm}"
    )

    # next step: replace this print-only stub with real pipeline calls
