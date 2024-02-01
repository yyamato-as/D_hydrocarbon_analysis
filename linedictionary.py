# all the data from CDMS
line_dict = {}
line_dict["MWC_480"] = {"CCD_fs1": {"spw": "0", 
                                    "qn": "4-3", 
                                    "restfreq": "288.4988349GHz", # the rest frequency of the group 1 is the strongest transition among hfs
                                    "hfs": {"hfs1": {"spw": "0", "qn": "4_55-3_45", "restfreq": "288.4947994GHz"}, # second and third quantum numbers are half integer
                                            "hfs2": {"spw": "0", "qn": "4_54-3_44", "restfreq": "288.4957734GHz"},
                                            "hfs3": {"spw": "0", "qn": "4_56-3_45", "restfreq": "288.4988349GHz"},
                                            "hfs4": {"spw": "0", "qn": "4_55-3_44", "restfreq": "288.4991049GHz"},
                                            "hfs5": {"spw": "0", "qn": "4_54-3_43", "restfreq": "288.4991082GHz"},
                                           },
                                   },
                        "CCD_fs2": {"spw": "0", 
                                    "qn": "4-3", 
                                    "restfreq": "288.5544885GHz", # the rest frequency of the group 2 is the strongest transition among hfs
                                    "hfs": {"hfs1": {"spw": "0", "qn": "4_45-3_34", "restfreq": "288.5544885GHz"}, # second and third quantum numbers are half integer
                                            "hfs2": {"spw": "0", "qn": "4_44-3_33", "restfreq": "288.5545644GHz"},
                                            "hfs3": {"spw": "0", "qn": "4_43-3_32", "restfreq": "288.5548170GHz"},
                                            "hfs4": {"spw": "0", "qn": "4_43-3_33", "restfreq": "288.5567527GHz"},
                                            "hfs5": {"spw": "0", "qn": "4_44-3_34", "restfreq": "288.5575810GHz"},
                                            },
                                   },
                        "DCO+": {"spw": "2", "qn": "4-3", "restfreq": "288.1438583GHz", "hfs": None},
                        "DCN":  {"spw": "6", "qn": "4-3", "restfreq": "289.6449170GHz", "hfs": None},
                        }