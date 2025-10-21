## New timing treatment in the `RealisticCaloDigi` processor (03.2021)

### Possible integration methods: `Standard` VS `ROC`

The RealisticCaloDigi processor includes two different timing treatment methods, also influencing the energy integration.

- The `Standard` integration sums up the energy of Monte Carlo contributions (MCCs) of *SimCalorimeterHits* within a certain time window (`t_min` < t < `t_max`). The time stamp of the reconstructed hit is set to the earliest MCC time of the simulated hit and the energy is set to the sum of the MCC energies within the time window.
- The `ROC` (ReadOut Chip) integration reflects a bit more how the electronics designed by the Omega group behaves. These chips include a fast shaper `t_fast` and a slow shaper `t_slow` value which are use in this algorithm to perform the energy integration and hit time estimate. Starting from the earliest MCC, the algorithm checks if the accumulated energy within `t_fast` passes the energy threshold. If not, it goes to the next MCC and performs the same check. If the threshold is reached, the hit time is set to the current MCC time and the energy is set to the accumulated MCC energy within `t_slow` from the current MCC. An absolute time cut is finally applied. If the hit time is in a certain time window (`t_min` < t < `t_max`), it is accepted, else discarded. If the hit is accepted, a random gaussian smearing is (optionally) applied to the hit time.

### Processor parameters

Here are listed the processor parameters acting on timing.

General parameters:

- `timingWindowMin`: The minimum hit time to accept a hit (t > `t_min`). Float, unit in `ns`, default is -10.
- `timingWindowMax`: The maximum hit time to accept a hit (t < `t_max`). Float, unit in `ns`, default is 100.
- `integrationMethod`: The energy integration and timing estimate method for digitization. Options are "Standard" (default) and "ROC".

Parameters for the "Standard" method **only**:

- `timingCorrectForPropagation`: Hit time correction `t_corr` for propagation computed as radial distance divided by `c`. Hit time becomes t - t_corr. Int, default is 0 (false).

Parameters for the "ROC" method **only**:

- `fastShaper`: The fast shaper integration time `t_fast`. Must be set in the steering file if the ROC method is chosen. Float, unit in `ns`, default is 0.
- `slowShaper`: The slow shaper integration time `t_slow`. Must be set in the steering file if the ROC method is chosen. Float, unit in `ns`, default is 0.
- `timingResolution`: Apply a random gaussian smearing `t_res` on the final hit time. The parameter value corresponds to the width of the gaussian function. Only applied if greater than 0. Float, unit in `ns`, default is 0 (no time smearing).

The default behavior of the digitizer is then:

- `integrationMethod`: "Standard"
- `timingWindowMin`: -10 ns
- `timingWindowMax`: 100 ns
- `timingCorrectForPropagation`: false


### Additional remarks on the ROC implementation

#### Pull request
The ROC implementation has been implementation through the `PR83` in the MarlinReco Github repository: https://github.com/iLCSoft/MarlinReco/pull/83. More detailed information, in particular the discussion with detector experts, can be found there.

#### Current ROC parameter values for known technologies (03.2021)

- AHCal:
  - `fastShaper`: 15 ns
  - `slowShaper`: 50 ns
  - `timingResolution`: 700 ps
- SiWEcal:
  - `fastShaper`: 90 ns. First estimate
  - `slowShaper`: 180 ns. First estimate
  - `timingResolution`: 700 ps. First guess from AHCAL value, as the TDC is similar. Not a real estimate





