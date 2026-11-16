import math
from podio import root_io

EXPECTED_HIT_AND_CLUSTER_WEIGHTS = {
    1: ([{"e⁻": 1.0}, {"e⁻": 1.0}], [{"e⁻": 1.0}]),
    2: ([{"γ": 1.0}], [{"γ": 1.0}]),
    3: ([{"backscattered n": 1.0}], [{"backscattered n": 1.0}]),
    4: (
        [{"e⁻": 1.0}, {"backscattered e⁻": 1.0}, {"backscattered e⁻": 1.0}],
        [{"e⁻": 1.0}, {"backscattered e⁻": 1.0}],
    ),
    5: ([{"μ⁺": 1.0}], [{"μ⁺": 1.0}]),
    6: ([{"π⁺": 1.0}, {"π⁺": 1.0}], [{"π⁺": 1.0}]),
    7: ([{"γ": 1.0}], [{"γ": 1.0}]),
    8: ([{"π⁻": 1.0}], [{"π⁻": 1.0}]),
    9: ([{"p": 1.0}], [{"p": 1.0}]),
    10: ([{"π⁺": 1.0}], [{"π⁺": 1.0}]),
    11: ([{"π⁺": 1.0}], [{"π⁺": 1.0}]),
    12: ([{"π⁺": 1.0}], [{"π⁺": 1.0}]),
    13: ([{"p": 1.0}], [{"p": 1.0}]),
    14: ([{"π⁺": 1.0}], [{"π⁺": 1.0}]),
    15: ([{"π⁺": 1.0}], [{"π⁺": 1.0}]),
    16: ([{"p": 1.0}], [{"p": 1.0}]),
    17: ([{"δ e⁻": 1.0}], [{"δ e⁻": 1.0}]),
    18: ([{"e⁻": 1.0}, {"e⁻": 1.0}], [{"e⁻": 1.0}]),
    19: ([{"e⁻": 0.25, "γ": 0.75}], [{"e⁻": 0.25, "γ": 0.75}]),
    20: ([{"π⁺": 1.0}, {"π⁺": 1.0}], [{"π⁺": 1.0}]),
    21: ([{"e⁻": 1.0}, {}], [{"e⁻": 0.6667}]),
    22: ([{"e⁻": 1.0}, {}], [{"e⁻": 0.6667}]),
    23: ([{"γ": 1.0}, {}, {}], [{"γ": 1.0}, {}]),
    24: ([{"n": 1.0}, {"n": 1.0}], [{"n": 1.0}]),
    25: ([{"n": 1.0}, {"γ": 1.0}], [{"n": 1.0}, {"γ": 1.0}]),
    26: ([{"e⁻": 0.75, "γ": 0.25}], [{"e⁻": 0.75, "γ": 0.25}]),
    27: ([{"π⁺": 1.0}], [{"π⁺": 1.0}]),
}

KNOWN_FAILING_CASES = []


def get_object_key(podio_object):
    object_id = podio_object.getObjectID()
    return (object_id.collectionID, object_id.index)


def sum_weights_per_source(links, sources, particle_names):
    source_keys = [get_object_key(source) for source in sources]
    weights = [{} for _ in sources]
    for link in links:
        source_key = get_object_key(link.getFrom())
        source_weights = weights[source_keys.index(source_key)]
        particle_name = particle_names[link.getTo().getObjectID().index]
        source_weights[particle_name] = (
            source_weights.get(particle_name, 0.0) + link.getWeight()
        )
    return weights


def get_mismatch(found_weights, expected_weights):
    if found_weights.keys() != expected_weights.keys():
        return "wrong particles"
    for name, weight in expected_weights.items():
        if not math.isclose(found_weights[name], weight, abs_tol=1e-4):
            return "wrong weights"
    return ""


reader = root_io.Reader("truth_links.edm4hep.root")
all_as_expected = True
for frame in reader.get("events"):
    case = frame.get_parameter("number")
    print(f"{case}. {frame.get_parameter('description')}")
    particle_names = frame.get_parameter("particle_names").split(";")
    hits = list(frame.get("ECalHits")) + list(frame.get("HCalHits"))
    clusters = list(frame.get("Clusters"))
    hit_weights = sum_weights_per_source(
        frame.get("CaloHitMCParticleLinks"), hits, particle_names
    )
    cluster_weights = sum_weights_per_source(
        frame.get("ClusterMCParticleLinks"), clusters, particle_names
    )
    expected_hit_weights, expected_cluster_weights = EXPECTED_HIT_AND_CLUSTER_WEIGHTS[
        case
    ]
    case_passes = True
    for kind, found, expected in [
        ("hit", hit_weights, expected_hit_weights),
        ("cluster", cluster_weights, expected_cluster_weights),
    ]:
        for number, (found_weights, expected_weights) in enumerate(
            zip(found, expected)
        ):
            mismatch = get_mismatch(found_weights, expected_weights)
            case_passes = case_passes and not mismatch
            marker = f"❌ {mismatch}" if mismatch else "✅"
            print(
                f"   {marker} {kind} {number}: found {found_weights}, expected {expected_weights}"
            )
    if case_passes == (case in KNOWN_FAILING_CASES):
        all_as_expected = False
        print(
            f"   case {case} {'passes but is' if case_passes else 'fails but is not'} in KNOWN_FAILING_CASES"
        )
exit(0 if all_as_expected else 1)
