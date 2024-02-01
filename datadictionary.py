data_path = "/works/yamato/"
data_dict = {}
data_dict["MWC_480"] = {
    "MAPS_Band3": [
        data_path
        + f"MAPS_data/MWC_480/Band3/self-calibrated/MWC_480_{i}_spectral_line.ms"
        for i in ["LB1", "LB2", "LB3", "LB4", "LB6"]
    ],
    "MAPS_Band3_wcont": [
        data_path
        + f"MAPS_data/MWC_480/Band3/self-calibrated/MWC_480_{i}_spectral_line_wcont.ms"
        for i in ["LB1", "LB2", "LB3", "LB4", "LB6"]
    ],
    "MAPS_Band6": [
        data_path
        + f"MAPS_data/MWC_480/Band6/self-calibrated/MWC_480_{i}_spectral_line.ms"
        for i in ["SB1", "SB2", "LB1", "LB2", "LB3", "LB4", "LB5", "LB6"]
    ],
    "MAPS_Band6_wcont": [
        data_path
        + f"MAPS_data/MWC_480/Band6/self-calibrated/MWC_480_{i}_spectral_line_wcont.ms"
        for i in ["SB1", "SB2", "LB1", "LB2", "LB3", "LB4", "LB5", "LB6"]
    ],
    "Band7": [
        f"/works/yamato/C2D_disk/v1_product/measurement_set/MWC_480_{i}_spectral_line.ms"
        for i in ["EB1", "EB2", "EB3", "EB4", "EB5"]
    ],
    "Band7_wcont": [
        f"/works/yamato/C2D_disk/v1_product/measurement_set/MWC_480_{i}_spectral_line_wcont.ms"
        for i in ["EB1", "EB2", "EB3", "EB4", "EB5"]
    ],
}
