# Modified from MAPS script for three objects targetted in this project

disk_dict = {}

# disk_dict['GM_Aur'] = {'name': 'GM Aur',
#                        'distance': 159, #source distance in pc, from Gaia Collaboration 2018,
#                        'incl': 53.21, #inclination in degrees, Huang et al. (2020)
#                        'PA': 57.17, #position angle in degrees, Huang et al. (2020)
#                        'PA_gofish': 57.17 + 0.0, #position angle in degrees, corrected for gofish
#                        'PA_diskprojection': 57.17 + 180.0, #position angle in degrees, corrected for diskprojection
#                        'Teff': 4350, #Kelvin, Macias et al. 2018
#                        'L_star': 1.2, #stellar luminosity in solar luminosities, Macias et al. 2018
#                        'M_star': 1.1, #stellar mass in solar masses, Macias et al. 2018
#                        'M_star_dynamical': 1.1,
#                        'logMdot': -8.1, #log of stellar accretion rate in solar masses/yr (GM Aur is variable!), Ingleby et al. 2013
#                        'v_sys': 5.61, #LSR systemic velocity [km/s]; Huang et al. (2020)
#                        '12CO_extent' : 14.0, # km/s, measured in v0 images
#                        'RA_center': '04h55m10.98834s', #fitted to continuum peak during self-cal
#                        'Dec_center': '+030.21.58.879285' #fitted to continuum peak during self-cal
# }

disk_dict['MWC_480'] = {'name': 'MWC 480',
                        'distance': 161.8, #[pc], from Gaia DR2
                        'incl': 37.0, #[degrees], Liu et al. (2019)
                        'PA': 148.0, #[degrees], Liu et al. (2019)
                        'PA_gofish': 148 + 180.0, #position angle in degrees, corrected for gofish
                        'PA_diskprojection': 148.0 + 0.0, #position angle in degrees, corrected for diskprojection
                        'Teff': 8250.0, #[K], Montesinos et al. (2009)
                        'L_star': 21.9, #[solar luminosity], Montesinos et al. (2009). L_star seems highly uncertain, error -7.8, +20.7 according to Montesinos+09. Liu+19 quote 17.3 LSun, Simon+01 quote 11.5 LSun
                        'M_star': 2.11, #[solar masses], Simon et al. (2019)
                        'M_star_dynamical': 2.11,
                        'logMdot': -6.9, #log of stellar accretion rate in solar masses/yr, Liu et al. (2019)
                        'v_sys': 5.1, # LSR systemic velocity in [km/s]; Pietu et al (2007)
                        '12CO_extent' : 16.0, # km/s, measured in v0 images
                        'RA_center' : '04h58m46.27938s',
                        'Dec_center' : '+029d50m36.399246s'
                       }

disk_dict['AS_209'] = {'name': 'AS 209',
                       'distance': 121, #source distance in pc, from Gaia Collaboration 2018,
                       'incl': 34.97, #inclination in degrees, Huang et al. 2018b
                       'PA': 85.76, #position angle in degrees, Huang et al. 2018b
                       'PA_gofish': 85.76 + 0.0, #position angle in degrees, corrected for gofish
                       'PA_diskprojection': 85.76 + 0.0, #position angle in degrees, corrected for diskprojection
                       'Teff': 4266, #Kelvin, Andrews et al. 2009
                       'L_star': 1.41, #stellar luminosity in solar luminosities, Andrews et al. 2009
                       'M_star': 0.83, #stellar mass in solar masses, Andrews et al. 2018b
                       'M_star_dynamical': 1.2, #13CO rotation profile investigated by Rich
                       'logMdot': -7.3, #log of stellar accretion rate in solar masses/yr, Salyk et al. 2013
                       'v_sys': 4.6, #LSR systemic velocity [km/s]; Huang et al. (2017)
                       '12CO_extent' : 17.0, # km/s, measured in v0 image
                       'RA_center' :'16h49m15.29404s' ,
                       'Dec_center' : '-014.22.09.092905'
}

# disk_dict['IM_Lup'] = {'name': "IM Lup",
#                        'distance': 158.0, #source distance in pc, from Gaia Collaboration 2018,
#                        'incl': 47.5, #inclination in degrees, Huang et al. 2018b
#                        'PA': 144.5, #position angle in degrees, Huang et al. 2018b
#                        'PA_gofish': 144.5 + 0.0, #position angle in degrees, corrected for gofish
#                        'PA_diskprojection': 144.5 + 0.0, #position angle in degrees, corrected for diskprojection
#                        'Teff': 4266.0, #Kelvin, Alcala et al. 2017
#                        'L_star': 2.57, #stellar luminosity in solar luminosities, Alcala et al. 2017
#                        'M_star': 0.89, #stellar mass in solar masses, Andrews et al. 2018b
#                        'M_star_dynamical': 1.1,
#                        'logMdot':-7.9, #log of stellar accretion rate in solar masses/yr,Alcala et al. 2017
#                        'v_sys': 4.5, #LSR systemic velocity [km/s]; Pinte et al. (2018)
#                        '12CO_extent' : 12.0, # km/s, measured in v0 image
#                        'RA_center' : '15h56m09.18683s',
#                        'Dec_center': '-37.56.06.579654'
# }

disk_dict['HD_163296'] = {'name': "HD 163296",
                          'distance': 101.0, #source distance in pc, from Gaia Collaboration 2018,
                          'incl': 46.7, #inclination in degrees, Huang et al. 2018b
                          'PA': 133.3, #position angle in degrees, Huang et al. 2018b
                          'PA_gofish': 133.3 + 180.0, #position angle in degrees, corrected for gofish
                          'PA_diskprojection': 133.3 + 180.0, #position angle in degrees, corrected for diskprojection
                          'Teff': 9332.0, #Kelvin, Fairlamb et al. 2015
                          'L_star': 17.0, #stellar luminosity in solar luminosities, Fairlamb et al. 2015
                          'M_star': 2.0, #stellar mass in solar masses, Andrews et al. 2018b
                          'M_star_dynamical':2.0,
                          'logMdot': -7.4, #log of stellar accretion rate in solar masses/yr, Fairlamb et al. 2015
                          'v_sys': 5.76, #LSR systemic velocity [km/s]; Teague et al. (2019)
                          '12CO_extent' : 21.0, # km/s, measured in v0 images
                          'RA_center' : '17h56m21.27703s',
                          'Dec_center' : '-021.57.22.640647'
}
