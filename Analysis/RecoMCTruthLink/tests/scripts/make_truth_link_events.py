import sys
import math
import podio
import edm4hep
from podio import root_io

CALIBRATION = 40.0
SYSTEMS = ["ECal", "HCal"]


def vector(position):
    return edm4hep.Vector3d(position[1], 0.0, position[0])


class Event:
    def __init__(self):
        self.mc_particles = edm4hep.MCParticleCollection()
        self.contributions = edm4hep.CaloHitContributionCollection()
        self.sim_hits = [edm4hep.SimCalorimeterHitCollection() for _ in SYSTEMS]
        self.calo_hits = [edm4hep.CalorimeterHitCollection() for _ in SYSTEMS]
        self.calo_to_sim_links = [
            edm4hep.CaloHitSimCaloHitLinkCollection() for _ in SYSTEMS
        ]
        self.clusters = edm4hep.ClusterCollection()
        self.particle_names = []
        self.next_cell_id = 0

    def particle(
        self,
        name,
        pdg,
        vertex,
        endpoint,
        parent=None,
        charge=0.0,
        generator_status=None,
        backscatter=False,
        decayed_in_tracker=False,
        decayed_in_calorimeter=False,
        vertex_is_not_endpoint_of_parent=False,
    ):
        particle = self.mc_particles.create()
        particle.setPDG(pdg)
        particle.setCharge(charge)
        default_generator_status = 1 if parent is None else 0
        particle.setGeneratorStatus(
            default_generator_status if generator_status is None else generator_status
        )
        vertex_vector = vector(vertex)
        endpoint_vector = vector(endpoint)
        particle.setVertex(vertex_vector)
        particle.setEndpoint(endpoint_vector)
        length = math.dist(vertex, endpoint) or 1.0
        direction = (
            (endpoint[0] - vertex[0]) / length,
            (endpoint[1] - vertex[1]) / length,
        )
        momentum = edm4hep.Vector3d(10.0 * direction[1], 0.0, 10.0 * direction[0])
        particle.setMomentum(momentum)
        particle.setBackscatter(backscatter)
        particle.setDecayedInTracker(decayed_in_tracker)
        particle.setDecayedInCalorimeter(decayed_in_calorimeter)
        particle.setVertexIsNotEndpointOfParent(vertex_is_not_endpoint_of_parent)
        if parent is not None:
            particle.addToParents(parent)
            parent.addToDaughters(particle)
        self.particle_names.append(name)
        return particle

    def new_cell_id(self):
        self.next_cell_id += 1
        return self.next_cell_id

    def sim_hit(self, contributions, system=0):
        sim_hit = self.sim_hits[system].create()
        sim_hit.setCellID(self.new_cell_id())
        for particle, energy in contributions:
            contribution = self.contributions.create()
            contribution.setPDG(particle.getPDG())
            contribution.setEnergy(energy)
            contribution.setParticle(particle)
            sim_hit.addToContributions(contribution)
        sim_hit.setEnergy(sum(energy for _, energy in contributions))
        return sim_hit

    def calo_hit(self, sim_hits, energy, system=0):
        calo_hit = self.calo_hits[system].create()
        calo_hit.setCellID(self.new_cell_id())
        calo_hit.setEnergy(energy)
        for sim_hit in sim_hits:
            link = self.calo_to_sim_links[system].create()
            link.setFrom(calo_hit)
            link.setTo(sim_hit)
            link.setWeight(1.0)
        return calo_hit

    def hit(self, contributions, system=0):
        sim_hit = self.sim_hit(contributions, system)
        return self.calo_hit([sim_hit], CALIBRATION * sim_hit.getEnergy(), system)

    def cluster(self, calo_hits, energy_correction=1.0, position=(1600, 0)):
        cluster = self.clusters.create()
        cluster.setEnergy(energy_correction * sum(hit.getEnergy() for hit in calo_hits))
        position_vector = edm4hep.Vector3f(position[1], 0.0, position[0])
        cluster.setPosition(position_vector)
        for calo_hit in calo_hits:
            cluster.addToHits(calo_hit)


def electron_from_ip(event):
    """An electron from the IP makes two hits."""
    electron = event.particle("e⁻", 11, (0, 0), (1800, 0), charge=-1.0)
    first_hit = event.hit([(electron, 0.01)])
    second_hit = event.hit([(electron, 0.02)])
    event.cluster([first_hit, second_hit])


def gun_photon_without_parents(event):
    """A particle-gun photon (generator status 0, no parents) makes a hit."""
    photon = event.particle("γ", 22, (0, 0), (1800, 0), generator_status=0)
    hit = event.hit([(photon, 0.01)])
    event.cluster([hit])


def backscattered_neutron(event):
    """A backscattered neutron makes a hit on return to the calorimeter."""
    neutron = event.particle("n", 2112, (0, 0), (1800, 0))
    backscattered_neutron = event.particle(
        "backscattered n",
        2112,
        (1600, 0),
        (1600, 400),
        parent=neutron,
        backscatter=True,
    )
    hit = event.hit([(backscattered_neutron, 0.01)])
    event.cluster([hit])


def photon_radiated_by_backscattered_electron(event):
    """A backscattered electron radiates a photon in the tracker; the photon makes two hits."""
    electron = event.particle("e⁻", 11, (0, 0), (1900, 0), charge=-1.0)
    backscattered_electron = event.particle(
        "backscattered e⁻",
        11,
        (1560, 20),
        (1600, 500),
        parent=electron,
        charge=-1.0,
        backscatter=True,
    )
    photon = event.particle(
        "γ",
        22,
        (1100, 300),
        (1800, 700),
        parent=backscattered_electron,
        vertex_is_not_endpoint_of_parent=True,
    )
    electron_hit = event.hit([(electron, 0.03)])
    first_photon_hit = event.hit([(photon, 0.01)])
    second_photon_hit = event.hit([(photon, 0.02)])
    event.cluster([electron_hit])
    event.cluster([first_photon_hit, second_photon_hit])


def muon_from_pion_decay_in_tracker(event):
    """A π⁺ decays in the tracker; the μ⁺ makes a hit."""
    pion = event.particle(
        "π⁺", 211, (0, 0), (1000, 0), charge=1.0, decayed_in_tracker=True
    )
    muon = event.particle("μ⁺", -13, (1000, 0), (1800, 100), parent=pion, charge=1.0)
    hit = event.hit([(muon, 0.01)])
    event.cluster([hit])


def secondary_of_pion_from_kaon_decay_in_tracker(event):
    """A K⁰ₛ decays in the tracker; the π⁺ interacts in the calorimeter and its kept secondary proton makes two
    hits."""
    kaon = event.particle("K⁰ₛ", 310, (0, 0), (700, 50), decayed_in_tracker=True)
    pion = event.particle("π⁺", 211, (700, 50), (1800, 300), parent=kaon, charge=1.0)
    event.particle("π⁻", -211, (700, 50), (1900, -150), parent=kaon, charge=-1.0)
    proton = event.particle(
        "p", 2212, (1800, 300), (2100, 450), parent=pion, charge=1.0
    )
    first_hit = event.hit([(proton, 0.01)])
    second_hit = event.hit([(proton, 0.02)])
    event.cluster([first_hit, second_hit])


def photon_converts_in_calorimeter(event):
    """A photon converts in the calorimeter; the conversion electron makes a hit."""
    photon = event.particle("γ", 22, (0, 0), (1600, 0))
    electron = event.particle(
        "e⁻", 11, (1600, 0), (1800, 30), parent=photon, charge=-1.0
    )
    hit = event.hit([(electron, 0.01)])
    event.cluster([hit])


def three_kept_shower_generations(event):
    """A π⁻ showers; a kept third-generation shower photon makes a hit."""
    pion = event.particle("π⁻", -211, (0, 0), (1600, 0), charge=-1.0)
    proton = event.particle("p", 2212, (1600, 0), (1700, 20), parent=pion, charge=1.0)
    neutron = event.particle("n", 2112, (1700, 20), (1800, -20), parent=proton)
    photon = event.particle("γ", 22, (1800, -20), (1900, 0), parent=neutron)
    hit = event.hit([(photon, 0.01)])
    event.cluster([hit])


def proton_knocked_out_by_pion(event, vertex, pion_endpoint=(1900, 0)):
    pion = event.particle("π⁺", 211, (0, 0), pion_endpoint, charge=1.0)
    proton_endpoint = (vertex[0] + 200, 50)
    proton = event.particle(
        "p",
        2212,
        vertex,
        proton_endpoint,
        parent=pion,
        charge=1.0,
        vertex_is_not_endpoint_of_parent=True,
    )
    hit = event.hit([(proton, 0.01)])
    event.cluster([hit])
    return pion


def sister_pi0_decayed_in_tracker(event):
    """A π⁺ knocks out a proton and continues, 70 mm before its endpoint; a π⁰ from the same vertex decayed in
    the tracker."""
    pion = proton_knocked_out_by_pion(event, (1480, 0), pion_endpoint=(1550, 0))
    event.particle(
        "π⁰", 111, (1480, 0), (1480, 0), parent=pion, decayed_in_tracker=True
    )


def sister_pi0_decayed_in_calorimeter(event):
    """A π⁺ knocks out a proton and continues, 300 mm before its endpoint; a π⁰ from the same vertex decayed in
    the calorimeter."""
    pion = proton_knocked_out_by_pion(event, (1600, 0))
    event.particle(
        "π⁰", 111, (1600, 0), (1600, 0), parent=pion, decayed_in_calorimeter=True
    )


def sister_backscatter(event):
    """A π⁺ knocks out a proton and continues, 300 mm before its endpoint; a neutron from the same vertex is
    backscattered."""
    pion = proton_knocked_out_by_pion(event, (1600, 0))
    event.particle("n", 2112, (1600, 0), (1400, -100), parent=pion, backscatter=True)


def pion_from_decay_in_calorimeter(event):
    """A K⁰ₗ decays in the calorimeter; the π⁺ knocks out a proton and continues, 300 mm before its endpoint."""
    kaon = event.particle(
        "K⁰ₗ", 130, (0, 0), (1550, 0), generator_status=2, decayed_in_calorimeter=True
    )
    pion = event.particle(
        "π⁺", 211, (1550, 0), (1900, 0), parent=kaon, charge=1.0, generator_status=1
    )
    proton = event.particle(
        "p",
        2212,
        (1600, 0),
        (1800, 50),
        parent=pion,
        charge=1.0,
        vertex_is_not_endpoint_of_parent=True,
    )
    hit = event.hit([(proton, 0.01)])
    event.cluster([hit])


def knock_out_800_mm_before_endpoint(event):
    """A π⁺ knocks out a proton and continues, 800 mm before its endpoint."""
    proton_knocked_out_by_pion(event, (1100, 0))


def knock_out_50_mm_before_endpoint(event):
    """A π⁺ knocks out a proton and continues, 50 mm before its endpoint."""
    proton_knocked_out_by_pion(event, (1850, 0))


def sister_pi0_0p05_mm_away(event):
    """A π⁺ knocks out a proton and continues, 300 mm before its endpoint; a π⁰ is produced 0.05 mm away."""
    pion = proton_knocked_out_by_pion(event, (1600, 0))
    event.particle(
        "π⁰", 111, (1600, 0.05), (1600, 0.05), parent=pion, decayed_in_calorimeter=True
    )


def sister_pi0_0p5_mm_away(event):
    """A π⁺ knocks out a proton and continues, 300 mm before its endpoint; a π⁰ is produced 0.5 mm away."""
    pion = proton_knocked_out_by_pion(event, (1600, 0))
    event.particle(
        "π⁰", 111, (1600, 0.5), (1600, 0.5), parent=pion, decayed_in_calorimeter=True
    )


def delta_electron_in_tracker(event):
    """A π⁺ interacts in the tracker; the proton emits a δ-electron in the tracker, which makes a hit."""
    pion = event.particle(
        "π⁺", 211, (0, 0), (800, 0), charge=1.0, decayed_in_tracker=True
    )
    proton = event.particle("p", 2212, (800, 0), (1800, 50), parent=pion, charge=1.0)
    delta_electron = event.particle(
        "δ e⁻",
        11,
        (1200, 20),
        (1700, 150),
        parent=proton,
        charge=-1.0,
        vertex_is_not_endpoint_of_parent=True,
    )
    hit = event.hit([(delta_electron, 0.01)])
    event.cluster([hit])


def one_sim_hit_two_calo_hits(event):
    """An electron makes one sim hit, which the digitiser shares 60/40 between two calo hits."""
    electron = event.particle("e⁻", 11, (0, 0), (1800, 0), charge=-1.0)
    sim_hit = event.sim_hit([(electron, 0.01)])
    first_calo_hit = event.calo_hit([sim_hit], 0.6 * CALIBRATION * 0.01)
    second_calo_hit = event.calo_hit([sim_hit], 0.4 * CALIBRATION * 0.01)
    event.cluster([first_calo_hit, second_calo_hit])


def two_sim_hits_one_calo_hit(event):
    """An electron (10 MeV) and a photon (30 MeV) make two sim hits, which the digitiser merges into one calo
    hit."""
    electron = event.particle("e⁻", 11, (0, -110), (1700, -110), charge=-1.0)
    photon = event.particle("γ", 22, (0, 110), (1700, 110))
    electron_sim_hit = event.sim_hit([(electron, 0.01)])
    photon_sim_hit = event.sim_hit([(photon, 0.03)])
    calo_hit = event.calo_hit([electron_sim_hit, photon_sim_hit], CALIBRATION * 0.04)
    event.cluster([calo_hit])


def pion_in_ecal_and_hcal(event):
    """A π⁺ makes one hit in the ECAL and one in the HCAL."""
    pion = event.particle("π⁺", 211, (0, 0), (2000, 0), charge=1.0)
    ecal_hit = event.hit([(pion, 0.01)], system=0)
    hcal_hit = event.hit([(pion, 0.02)], system=1)
    event.cluster([ecal_hit, hcal_hit])


def noise_hit_in_corrected_cluster(event):
    """An electron cluster of a 0.4 GeV hit and a 0.2 GeV noise hit, with cluster energy corrected ×1.3."""
    electron = event.particle("e⁻", 11, (0, 0), (1800, 0), charge=-1.0)
    electron_hit = event.hit([(electron, 0.01)])
    noise_hit = event.calo_hit([], 0.2)
    event.cluster([electron_hit, noise_hit], energy_correction=1.3)


def noise_hit_in_uncorrected_cluster(event):
    """An electron cluster of a 0.4 GeV hit and a 0.2 GeV noise hit."""
    electron = event.particle("e⁻", 11, (0, 0), (1800, 0), charge=-1.0)
    electron_hit = event.hit([(electron, 0.01)])
    noise_hit = event.calo_hit([], 0.2)
    event.cluster([electron_hit, noise_hit])


def noise_only_cluster(event):
    """A cluster of noise hits only, 0.1 rad in polar angle from a photon."""
    photon = event.particle("γ", 22, (0, 0), (1800, 150))
    photon_hit = event.hit([(photon, 0.01)])
    first_noise_hit = event.calo_hit([], 0.2)
    second_noise_hit = event.calo_hit([], 0.1)
    event.cluster([photon_hit], position=(1650, 137))
    event.cluster([first_noise_hit, second_noise_hit], position=(1650, 300))


def backscattered_neutron_captured_in_tracker(event):
    neutron = event.particle("n", 2112, (0, 0), (1900, 0))
    backscattered_neutron = event.particle(
        "backscattered n",
        2112,
        (1700, 0),
        (1200, 200),
        parent=neutron,
        backscatter=True,
        decayed_in_tracker=True,
    )
    photon = event.particle(
        "γ", 22, (1200, 200), (1700, 150), parent=backscattered_neutron
    )
    neutron_hit = event.hit([(neutron, 0.03)])
    photon_hit = event.hit([(photon, 0.01)])
    return neutron_hit, photon_hit


def capture_photon_in_neutron_cluster(event):
    """A backscattered neutron is captured in the tracker; the capture photon's hit joins the neutron's cluster."""
    neutron_hit, photon_hit = backscattered_neutron_captured_in_tracker(event)
    event.cluster([neutron_hit, photon_hit])


def capture_photon_in_own_cluster(event):
    """A backscattered neutron is captured in the tracker; the capture photon's hit forms its own cluster."""
    neutron_hit, photon_hit = backscattered_neutron_captured_in_tracker(event)
    event.cluster([neutron_hit])
    event.cluster([photon_hit])


def two_particles_in_one_sim_hit(event):
    """An electron (10 + 20 MeV) and a photon (10 MeV) contribute to one sim hit."""
    electron = event.particle("e⁻", 11, (0, -50), (1800, -50), charge=-1.0)
    photon = event.particle("γ", 22, (0, 50), (1800, 50))
    hit = event.hit([(electron, 0.01), (electron, 0.02), (photon, 0.01)])
    event.cluster([hit])


def pion_and_its_shower_photon_in_one_sim_hit(event):
    """A π⁺ and its own shower photon contribute to one sim hit."""
    pion = event.particle("π⁺", 211, (0, 0), (1800, 0), charge=1.0)
    photon = event.particle("γ", 22, (1600, 0), (1700, 20), parent=pion)
    hit = event.hit([(pion, 0.01), (photon, 0.01)])
    event.cluster([hit])


EVENTS = [
    electron_from_ip,
    gun_photon_without_parents,
    backscattered_neutron,
    photon_radiated_by_backscattered_electron,
    muon_from_pion_decay_in_tracker,
    secondary_of_pion_from_kaon_decay_in_tracker,
    photon_converts_in_calorimeter,
    three_kept_shower_generations,
    sister_pi0_decayed_in_tracker,
    sister_pi0_decayed_in_calorimeter,
    sister_backscatter,
    pion_from_decay_in_calorimeter,
    knock_out_800_mm_before_endpoint,
    knock_out_50_mm_before_endpoint,
    sister_pi0_0p05_mm_away,
    sister_pi0_0p5_mm_away,
    delta_electron_in_tracker,
    one_sim_hit_two_calo_hits,
    two_sim_hits_one_calo_hit,
    pion_in_ecal_and_hcal,
    noise_hit_in_corrected_cluster,
    noise_hit_in_uncorrected_cluster,
    noise_only_cluster,
    capture_photon_in_neutron_cluster,
    capture_photon_in_own_cluster,
    two_particles_in_one_sim_hit,
    pion_and_its_shower_photon_in_one_sim_hit,
]

writer = root_io.Writer(sys.argv[1])
for number, build_event in enumerate(EVENTS, start=1):
    event = Event()
    build_event(event)
    description = " ".join(build_event.__doc__.split())
    particle_names = ";".join(event.particle_names)
    frame = podio.Frame()
    frame.put_parameter("number", number)
    frame.put_parameter("name", build_event.__name__)
    frame.put_parameter("description", description)
    frame.put_parameter("particle_names", particle_names)
    frame.put(event.mc_particles, "MCParticles")
    frame.put(event.contributions, "CaloHitContributions")
    for system, name in enumerate(SYSTEMS):
        frame.put(event.sim_hits[system], f"{name}SimHits")
        frame.put(event.calo_hits[system], f"{name}Hits")
        frame.put(event.calo_to_sim_links[system], f"{name}Links")
    frame.put(event.clusters, "Clusters")
    frame.put(edm4hep.TrackCollection(), "EmptyTracks")
    writer.write_frame(frame, "events")
del writer
