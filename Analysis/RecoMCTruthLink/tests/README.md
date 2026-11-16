# RecoMCTruthLinker tests

27 synthetic events, one per case, test the calo hit → MCParticle and cluster → MCParticle links. The same events and check run for `CreateTruthLinks` in k4RecCalorimeter; both linkers must give identical links.

`make_truth_link_events.py` writes the events → `runRecoMCTruthLinker.py` runs the linker with `k4run` → `check_truth_links.py` compares every link with the expected values below. Cases in `KNOWN_FAILING_CASES` must fail, all others must pass.

## Weights

- Calo hit: fraction of the calo hit energy from the MCParticle. A calo hit made of several sim hits is shared among them by sim hit energy.
- Cluster: calo hit energy from the MCParticle over the summed calo hit energy of the cluster, noise hits included.

## Originator

Each hit is linked to the particle that entered the calorimeter (the originator), found by walking up the parents of the particle in the sim hit contribution:

- A generator particle, or one without parents, is its own originator.
- A backscattered particle is its own originator.
- A particle whose parent decayed in the tracker is its own originator.
- A particle from a non-destructive interaction (parent continues) is its own originator if the interaction was in the tracker. Decided by a sister at the same vertex (within 0.1 mm): decayed in the tracker → tracker; π⁰ or backscatter → calorimeter. Otherwise by the grandparent decayed in the calorimeter → calorimeter, else by the distance from the vertex to the parent's endpoint: > 200 mm → tracker.
- Otherwise the parent's originator.

In the cluster loop, the hits of a particle with a backscattered parent are moved to an ancestor in the same cluster, if there is one.

## Cases

| # | Event | Calo hit links | Cluster links | Why |
|---|---|---|---|---|
| 1 | An electron from the IP makes two hits. | e⁻ 1; e⁻ 1 | e⁻ 1 | generator particle |
| 2 | A particle-gun photon (generator status 0, no parents) makes a hit. | γ 1 | γ 1 | no parents |
| 3 | A backscattered neutron makes a hit on return to the calorimeter. | backscattered n 1 | backscattered n 1 | backscatter |
| 4 | A backscattered electron radiates a photon in the tracker; the photon makes two hits. | e⁻ 1; backscattered e⁻ 1; backscattered e⁻ 1 | e⁻ 1; backscattered e⁻ 1 | photon → its backscattered ancestor |
| 5 | A π⁺ decays in the tracker; the μ⁺ makes a hit. | μ⁺ 1 | μ⁺ 1 | parent decayed in the tracker |
| 6 | A K⁰ₛ decays in the tracker; the π⁺ interacts in the calorimeter and its secondary proton makes two hits. | π⁺ 1; π⁺ 1 | π⁺ 1 | π⁺ entered the calorimeter |
| 7 | A photon converts in the calorimeter; the conversion electron makes a hit. | γ 1 | γ 1 | conversion in the calorimeter |
| 8 | A π⁻ showers; a third-generation shower photon makes a hit. | π⁻ 1 | π⁻ 1 | whole shower in the calorimeter |
| 9 | A π⁺ knocks out a proton and continues, 70 mm before its endpoint; a π⁰ from the same vertex decayed in the tracker. | p 1 | p 1 | sister decayed in the tracker |
| 10 | As 9, 300 mm before the endpoint; the π⁰ decayed in the calorimeter. | π⁺ 1 | π⁺ 1 | π⁰ sister |
| 11 | As 9, 300 mm before the endpoint; a backscattered neutron from the same vertex. | π⁺ 1 | π⁺ 1 | backscattered sister |
| 12 | A K⁰ₗ decays in the calorimeter; the π⁺ knocks out a proton and continues, 300 mm before its endpoint. | π⁺ 1 | π⁺ 1 | grandparent decayed in the calorimeter |
| 13 | A π⁺ knocks out a proton and continues, 800 mm before its endpoint. | p 1 | p 1 | > 200 mm |
| 14 | As 13, 50 mm before the endpoint. | π⁺ 1 | π⁺ 1 | < 200 mm |
| 15 | As 13, 300 mm before the endpoint; a π⁰ is produced 0.05 mm away. | π⁺ 1 | π⁺ 1 | π⁰ sister at the same vertex |
| 16 | As 15, the π⁰ 0.5 mm away. | p 1 | p 1 | π⁰ not at the same vertex; > 200 mm |
| 17 | A π⁺ interacts in the tracker; the proton emits a δ-electron in the tracker, which makes a hit. | δ e⁻ 1 | δ e⁻ 1 | > 200 mm |
| 18 | An electron makes one sim hit, which the digitiser shares 60/40 between two calo hits. | e⁻ 1; e⁻ 1 | e⁻ 1 | both calo hits linked |
| 19 | An electron (10 MeV) and a photon (30 MeV) make two sim hits, which the digitiser merges into one calo hit. | e⁻ 0.25, γ 0.75 | e⁻ 0.25, γ 0.75 | shared by sim hit energy |
| 20 | A π⁺ makes one hit in the ECAL and one in the HCAL. | π⁺ 1; π⁺ 1 | π⁺ 1 | both calorimeters |
| 21 | An electron cluster of a 0.4 GeV hit and a 0.2 GeV noise hit, with cluster energy corrected ×1.3. | e⁻ 1; none | e⁻ 0.67 | normalised to summed calo hit energy, not cluster energy |
| 22 | As 21, without correction. | e⁻ 1; none | e⁻ 0.67 | noise hit counts in the denominator |
| 23 | A cluster of noise hits only, 0.1 rad in polar angle from a photon; the photon makes its own cluster. | γ 1; none; none | γ 1; none | noise hits are not linked |
| 24 | A backscattered neutron is captured in the tracker; the capture photon's hit joins the neutron's cluster. | n 1; n 1 | n 1 | photon moved to its ancestor in the same cluster |
| 25 | As 24, the capture photon's hit forms its own cluster. | n 1; γ 1 | n 1; γ 1 | no ancestor in the photon's cluster |
| 26 | An electron (10 + 20 MeV) and a photon (10 MeV) contribute to one sim hit. | e⁻ 0.75, γ 0.25 | e⁻ 0.75, γ 0.25 | shared by contribution energy |
| 27 | A π⁺ and its own shower photon contribute to one sim hit. | π⁺ 1 | π⁺ 1 | both contributions → π⁺ |

"; " separates calo hits or clusters, in the order of the event; "none" means no link.
