import math
import pathlib
import numpy as np
import sinter
import stim
from ldpc import bp_decoder

class BeliefPropagation(sinter.Decoder):
    def decode_via_files(self,
                         *,
                         num_shots: int,
                         num_dets: int,
                         num_obs: int,
                         dem_path: pathlib.Path,
                         dets_b8_in_path: pathlib.Path,
                         obs_predictions_b8_out_path: pathlib.Path,
                         tmp_dir: pathlib.Path,
                       ) -> None:
        # Read input data. Note it's not a great idea to store it in memory all at once like this.
        detector_error_model = stim.DetectorErrorModel.from_file(dem_path)
        packed_detection_event_data = np.fromfile(dets_b8_in_path, dtype=np.uint8)
        packed_detection_event_data.shape = (num_shots, math.ceil(num_dets / 8))

        # Make predictions.
        all_predictions = []
        for shot in packed_detection_event_data:
            unpacked = np.unpackbits(shot, bitorder='little')
            shot_predictions = []

            # This is a terrible way to make the prediction.
            # It doesn't even look at the error model!
            excitations = np.count_nonzero(unpacked)
            for k in range(num_obs):
                shot_predictions.append(excitations % (k + 2) == 0)

            all_predictions.append(shot_predictions)

        # Write predictions.
        np.packbits(all_predictions, axis=1, bitorder='little').tofile(obs_predictions_b8_out_path)

    def compile_decoder_for_dem(self, *, dem: stim.DetectorErrorModel) -> 'sinter.CompiledDecoder':
        # This will be added in v1.12 and sinter will prefer it, as it avoids the disk as a bottleneck.
        #
        # You have to return an object with this method:
        #    def decode_shots_bit_packed(
        #                self,
        #                *,
        #                bit_packed_detection_event_data: np.ndarray,
        #        ) -> np.ndarray:
        raise NotImplementedError()