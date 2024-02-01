# all the data from CDMS
line_dict = {}
line_dict["MWC_480"] = {
    "C2H_1-0_fs1": {
        "data": "MAPS_Band3",
        "Band": "3",
        "spw": "0",
        "qn": "1_2-0_1",
        "restfreq": "87.3168980GHz",
        "freqs": ["87.284108GHz", "87.316894GHz", "87.328587GHz"],
    },
    # "C2H_1-0_fs2": {
    #     "spw": "0",
    #     "qn": "1_21-0_10",
    #     "restfreq": "87.3285850GHz",
    # },
    "C2H_3-2": {
        "data": "MAPS_Band6",
        "Band": "6",
        "spw": "3",
        "qn": "3-2",
        "restfreq": "262.0042600GHz",  # restfreq of the brightest component; see Guzman+21
        "freqs": [
            "262.0042600GHz",
            "262.0064820GHz",
            "262.0649860GHz",
            "262.0674690GHz",
            "262.0789347GHz",
        ],
    },
    "C2H_3-2_fs1": {
        "data": "MAPS_Band6",
        "Band": "6",
        "spw": "3",
        "qn": "3_4-2_3",
        "restfreq": "262.0042600GHz",
        "freqs": [
            "262.0042600GHz",
            "262.0064820GHz",
        ],
    },
    "C2H_3-2_fs2": {
        "data": "MAPS_Band6",
        "Band": "6",
        "spw": "3",
        "qn": "3_3-2_2",
        "restfreq": "262.0649860GHz",
        "freqs": [
            "262.0649860GHz",
            "262.0674690GHz",
            "262.0789347GHz",
        ],
    },
    "C2D_4-3": {
        "data": "Band7",
        "Band": "7",
        "spw": "0",
        "qn": "4-3",
        "restfreq": "288.4988349GHz",  # the rest frequency of the group 1 is the strongest transition among hfs
        "freqs": [
            "288.4988349GHz",
            "288.4991049GHz",
            "288.4991082GHz",
            "288.5544885GHz",
            "288.5545644GHz",
            "288.5548170GHz",
        ],
    },
    "C2D_4-3_fs1": {
        "data": "Band7",
        "Band": "7",
        "spw": "0",
        "qn": "4_5-3_4",
        "restfreq": "288.4988349GHz",  # the rest frequency of the group 1 is the strongest transition among hfs
        "freqs": ["288.4988349GHz", "288.4991049GHz", "288.4991082GHz"],
    },
    "C2D_4-3_fs2": {
        "data": "Band7",
        "Band": "7",
        "spw": "0",
        "qn": "4_4-3_3",
        "restfreq": "288.5544885GHz",  # the rest frequency of the group 2 is the strongest transition among hfs
        "freqs": ["288.5544885GHz", "288.5545644GHz", "288.5548170GHz"],
    },
    "c-C3HD_9-8": {
        "spw": "1",
        "qn": "9-8",
        "restfreq": "286.9333850GHz",  # the rest frequency of the group 2 is the strongest transition among hfs
        "maskfreqs": [
            "286.9333465GHz",
            "286.9333850GHz",
            "286.9334115GHz",
            "286.9333850GHz",
        ],
    },
    "DCO+_4-3": {
        "spw": "2",
        "qn": "4-3",
        "restfreq": "288.1438583GHz",
        "maskfreqs": ["288.1438583GHz"],
    },
    "H2CO_4_04-3_03": {
        "spw": "3",
        "qn": "4_04-3_03",
        "restfreq": "290.6234050GHz",
        "maskfreqs": ["290.6234050GHz"],
    },
    "C34S_6-5": {
        "spw": "4",
        "qn": "6-5",
        "restfreq": "289.2090684GHz",
        "maskfreqs": ["289.2090684GHz"],
    },
    "CH3OH_6_06-5_05": {
        "spw": "4",
        "qn": "6-5",
        "restfreq": "289.2090684GHz",
        "maskfreqs": ["289.2090684GHz"],
    },
    "DCN_4-3": {
        "spw": "6",
        "qn": "4-3",
        "restfreq": "289.6449170GHz",
        "maskfreqs": ["289.6449170GHz"],
    },
    "SO_7-6": {
        "spw": "8",
        "qn": "3Sigma 7_7-6_6",
        "restfreq": "301.2861240GHz",
        "maskfreqs": ["301.2861240GHz"],
    },
    "H2CO_4_13-3_12": {
        "spw": "9",
        "qn": "4_13-3_12",
        "restfreq": "300.8366350GHz",
        "maskfreqs": ["300.8366350GHz"],
    },
}
