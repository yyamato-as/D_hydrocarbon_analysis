import xmltodict
import xml.etree.ElementTree as ET
import astropy.constants as ac
import numpy as np
import os
from scipy.interpolate import interp1d

h = ac.h.cgs.value
c = ac.c.cgs.value
k_B = ac.k_B.cgs.value


class QuantumNumber:
    def __init__(self, **qns):
        for qn in qns:
            setattr(self, qn, qns[qn])


class State:
    def __init__(self, energy, weight, QN=None):
        """_summary_

        Parameters
        ----------
        energy : float
            energy of the level in cm-1
        weight : float
            total statistical weight of the level
        QN : QuantumNumber class, optional
            quantum number, by default None
        """
        self.energy = energy
        self.weight = weight
        self.QN = QN


class RadiativeTransition:
    def __init__(
        self,
        species,
        restfreq,
        restfreq_err,
        gup,
        glow,
        EinsteinA,
        Eup=None,
        Elow=None,
        QNup=None,
        QNlow=None,
        comment=None,
    ):
        self.species = species
        self.restfreq = restfreq  # in Hz
        self.restfreq_err = restfreq_err  # in Hz
        self.gup = gup
        self.glow = glow
        self.EinsteinA = EinsteinA

        # energy levels
        if Eup is None:
            if Elow is None:
                raise ValueError("Please provide either Eup or Elow.")
            self.Elow = Elow
            self.Eup = self.Elow + h * self.restfreq / k_B
        else:
            self.Eup = Eup
            self.Elow = self.Eup - h * self.restfreq / k_B

        self.QNup = QNup
        self.QNlow = QNlow
        self.comment = comment


class CollisionalTransition:
    def __init__(
        self,
        species,
        temperature,
        coefficient,
        collider=None,
        QNup=None,
        QNlow=None,
        comment=None,
    ):
        self.species = species
        self.collider = collider
        self.temperature = temperature
        self.coefficient = coefficient
        self.QNup = QNup
        self.QNlow = QNlow
        self.comment = comment


class Transitions:
    def __init__(
        self, transitions, states=None, lowerid=None, upperid=None, colpartners=None
    ):
        self.states = states
        self.transitions = transitions
        self.lowerid = lowerid
        self.upperid = upperid
        self.colpartners = colpartners

        self.nstates = len(states)
        self.ntrans = len(transitions)


class PartitionFunction:
    def __init__(self, specie, T, Q, database=None, ntrans=None):
        self.specie = specie
        self.T = T
        self.Q = Q
        self.database = database
        self.ntrans = ntrans

        self.function = self._get_function()

    def __call__(self, T):
        val = self.function(T)
        if val.size == 1:
            return float(val)
        else:
            return val

    def _get_function(self):
        return interp1d(self.T, self.Q, kind="cubic", fill_value="extrapolate")


def wavenumber_to_Kelvin(wavenumber):
    return wavenumber * h * c / k_B


def Kelvin_to_wavenumber(Kelvin):
    return Kelvin * k_B / (h * c)


class XSAMSFile:
    def __init__(self, filename):
        self.filename = filename
        xmldata = ET.parse(self.filename).getroot()
        self.datadict = xmltodict.parse(ET.tostring(xmldata))["ns0:XSAMSData"]
        self.molecules = self.datadict["ns0:Species"]["ns0:Molecules"]["ns0:Molecule"]
        if not isinstance(self.molecules, list):
            self.molecules = [self.molecules]

    def get_species_dict(self):
        IDs = [d["@speciesID"] for d in self.molecules]
        comments = [
            d["ns0:MolecularChemicalSpecies"]["ns0:Comment"] for d in self.molecules
        ]
        return {id: com for id, com in zip(IDs, comments)}

    def select_species(self, speciesID=None):
        self.speciesID = speciesID
        if self.speciesID is None:
            if len(self.molecules) > 1:
                raise ValueError(
                    "Provide speciesID to select one from multiple entries."
                )
            else:
                self.speciesID = self.molecules[0]["@speciesID"]
        ids = [d["@speciesID"] for d in self.molecules]
        self.molecule = self.molecules[ids.index(self.speciesID)]

    def load_partition_function(self):
        pf = self.molecule["ns0:MolecularChemicalSpecies"]["ns0:PartitionFunction"]
        self.T = list(map(float, pf["ns0:T"]["ns0:DataList"].split()))
        self.Q = list(map(float, pf["ns0:Q"]["ns0:DataList"].split()))

    def save_partition_function(self, filename, **kwargs):
        np.savetxt(filename, np.vstack([self.T, self.Q]).T, header="T(K)\t Q", **kwargs)

    def load_molecular_properties(self):
        # species ID used VAMDC-wide
        self.VAMDCSpeciesID = self.molecule["ns0:MolecularChemicalSpecies"][
            "ns0:VAMDCSpeciesID"
        ]

        # species name
        self.species = self.molecule["ns0:MolecularChemicalSpecies"][
            "ns0:ChemicalName"
        ]["ns0:Value"]

        # molecular formula
        self.formula = self.molecule["ns0:MolecularChemicalSpecies"][
            "ns0:OrdinaryStructuralFormula"
        ]["ns0:Value"]

        # molecular weight
        try:
            self.molecular_weight = float(
                self.molecule["ns0:MolecularChemicalSpecies"][
                    "ns0:StableMolecularProperties"
                ]["ns0:MolecularWeight"]["ns0:Value"]["#text"]
            )
        except:
            self.molecular_weight = None

    def load_state_data(self):
        # states
        self.stateid = []
        self.states = []
        for state in self.molecule["ns0:MolecularState"]:
            # skip the energy origin that is doubly entried
            if "@auxillary" in state and state["@auxillary"] == "true":
                continue
            # stateid
            self.stateid.append(state["@stateID"])

            # quantum numbers
            qn_dict = {}
            key = [key for key in state["ns0:Case"].keys() if "QNs" in key][0]
            for qn, val in state["ns0:Case"][key].items():
                qn_dict[qn.split(":")[1]] = (
                    val["#text"] if isinstance(val, dict) else val
                )
            QN = QuantumNumber(**qn_dict)

            # energy and weight
            property = state["ns0:MolecularStateCharacterisation"]
            energy = float(property["ns0:StateEnergy"]["ns0:Value"]["#text"])

            try:
                weight = float(property["ns0:TotalStatisticalWeight"])
            except KeyError:
                weight = None

            # state
            s = State(energy=energy, weight=weight, QN=QN)
            self.states.append(s)

        # number of energy levels
        self.nstates = len(self.states)

        # print(self.stateid)

        # if "ns0:Radiative" in self.datadict["ns0:Processes"]:
        #     self.load_radiative_transitions()
        # elif "ns0:Collisions" in self.datadict["ns0:Processes"]:
        #     self.load_collisional_transitions()
        # else:
        #     pass

    def load_radiative_transitions(self):
        # transition
        self.upperid = []
        self.lowerid = []
        self.transitions = []
        for trans in self.datadict["ns0:Processes"]["ns0:Radiative"][
            "ns0:RadiativeTransition"
        ]:
            if self.speciesID is not None and trans["ns0:SpeciesRef"] != self.speciesID:
                continue
            # restfreq
            restfreqs = trans["ns0:EnergyWavelength"]["ns0:Frequency"]
            if not isinstance(restfreqs, list):
                restfreqs = [restfreqs]
            accuracy = np.inf
            for nu0 in restfreqs:
                if float(nu0["ns0:Accuracy"]) >= accuracy:
                    continue
                restfreq = float(nu0["ns0:Value"]["#text"])
                restfreq_err = float(nu0["ns0:Accuracy"])

            # lower and upper state
            self.lowerid.append(trans["ns0:LowerStateRef"])
            self.upperid.append(trans["ns0:UpperStateRef"])

            # Einstein A
            EinsteinA = float(
                trans["ns0:Probability"]["ns0:TransitionProbabilityA"]["ns0:Value"][
                    "#text"
                ]
            )

            # state data
            upperstate = self.states[self.stateid.index(trans["ns0:UpperStateRef"])]
            lowerstate = self.states[self.stateid.index(trans["ns0:LowerStateRef"])]
            gup = upperstate.weight
            Eup = wavenumber_to_Kelvin(upperstate.energy)
            QNup = upperstate.QN
            glow = lowerstate.weight
            Elow = wavenumber_to_Kelvin(lowerstate.energy)
            QNlow = lowerstate.QN

            t = RadiativeTransition(
                species=self.species,
                restfreq=restfreq,
                restfreq_err=restfreq_err,
                gup=gup,
                glow=glow,
                EinsteinA=EinsteinA,
                Eup=Eup,
                Elow=Elow,
                QNup=QNup,
                QNlow=QNlow,
            )
            self.transitions.append(t)

        # number of transitions
        self.ntrans = len(self.transitions)

    def load_collisional_transitions(self, colpartnerIDs=None, sourceID=None):
        # transition
        self.colpartnerIDs = (
            colpartnerIDs  # for collider selection for collisional transitions
        )
        self.sourceID = sourceID  # source of collisional coefficients
        self.colpartners = []
        self.upperid = []
        self.lowerid = []
        self.transitions = []

        for trans in self.datadict["ns0:Processes"]["ns0:Collisions"][
            "ns0:CollisionalTransition"
        ]:
            # check source
            sources = trans["ns0:SourceRef"]
            if self.sourceID is not None and not self.sourceID in sources:
                continue

            # get colliders and upper and lower states
            for r, p in zip(trans["ns0:Reactant"], trans["ns0:Product"]):
                if r == p:
                    colpartID = r["ns0:SpeciesRef"]
                else:
                    upperstate = r
                    lowerstate = p

            if self.colpartnerIDs is not None and not colpartID in self.colpartnerIDs:
                continue

            if self.speciesID is not None and (
                upperstate["ns0:SpeciesRef"] != self.speciesID
                or lowerstate["ns0:SpeciesRef"] != self.speciesID
            ):
                continue

            self.colpartners.append(colpartID)

            # upper and lower states
            self.upperid.append(upperstate["ns0:StateRef"])
            self.lowerid.append(lowerstate["ns0:StateRef"])

            upperstate = self.states[self.stateid.index(upperstate["ns0:StateRef"])]
            lowerstate = self.states[self.stateid.index(lowerstate["ns0:StateRef"])]
            QNup = upperstate.QN
            QNlow = lowerstate.QN

            # temperature and coefficient
            data = trans["ns0:DataSets"]["ns0:DataSet"]["ns0:TabulatedData"]
            T = list(map(float, data["ns0:X"]["ns0:DataList"].split()))
            coeff = list(map(float, data["ns0:Y"]["ns0:DataList"].split()))
            comment = trans["ns0:Comments"]

            t = CollisionalTransition(
                species=self.species,
                temperature=T,
                coefficient=coeff,
                QNup=QNup,
                QNlow=QNlow,
                comment=comment,
            )

            self.transitions.append(t)

        # number of collision partners
        self.ncolpartners = len(set(self.colpartners))
        # number of transitions
        self.ntrans = len(self.transitions)


class SpectroscopicData:
    def __init__(self):
        pass

    def read_radiative_XSAMS_file(self, filename, speciesID=None):
        self.radiative = XSAMSFile(filename=filename)
        self.radiative.select_species(speciesID=speciesID)
        self.radiative.load_molecular_properties()
        self.radiative.load_state_data()
        self.radiative.load_radiative_transitions()

    def read_collisional_XSAMS_file(
        self, filename, speciesID=None, colpartnerIDs=None, sourceID=None
    ):
        self.collisional = XSAMSFile(
            filename=filename,
        )
        self.collisional.select_species(speciesID=speciesID)
        self.collisional.load_molecular_properties()
        self.collisional.load_state_data()
        self.collisional.load_collisional_transitions(
            colpartnerIDs=colpartnerIDs, sourceID=sourceID
        )

    def crossmatch_data(self):
        # first assert the data is for the same molecules
        assert (
            self.radiative.VAMDCSpeciesID == self.collisional.VAMDCSpeciesID
        ), "The molecular species do not match between radiative and collisional data"

        # intergate molecular properties
        self.species = self.radiative.species
        self.formula = self.radiative.formula
        self.molecular_weight = self.collisional.molecular_weight

        # assume that all of the states from collsional data are included in the radiative data
        ref_QNs = [s.QN for s in self.collisional.states]

        # modify the radiative data
        rad_stateid_new = []
        rad_states_new = []
        for id, s in zip(self.radiative.stateid, self.radiative.states):
            # # first check the energy difference to get the candidate
            # ref_idx = np.argmin(np.abs(np.array(ref_energies) - s.energy))
            # ref_idx = int(ref_idx)
            # print(s.energy, ref_energies[ref_idx])
            # if not np.isclose(s.energy, ref_energies[ref_idx], rtol=1e-3, atol=1e-8):
            #     matched = False
            # else:
            # qunatum number matching
            for i, QN in enumerate(ref_QNs):
                matched = True
                for qn, val in s.QN.__dict__.items():
                    if hasattr(QN, qn):
                        ref_val = getattr(QN, qn)
                        try:
                            val = float(val)
                            ref_val = float(ref_val)
                        except ValueError:
                            pass
                        if val != ref_val:
                            matched = False
                if matched:
                    ref_idx = i
                    break

            if matched:
                # modify radiative data
                rad_stateid_new.append(id)
                rad_states_new.append(s)
                self.collisional.upperid = [
                    id if uid == self.collisional.stateid[ref_idx] else uid
                    for uid in self.collisional.upperid
                ]
                self.collisional.lowerid = [
                    id if lid == self.collisional.stateid[ref_idx] else lid
                    for lid in self.collisional.lowerid
                ]
                self.collisional.stateid[ref_idx] = id
                self.collisional.states[ref_idx] = s
            else:
                # remove radiative data
                self.radiative.transitions = [
                    trans
                    for i, trans in enumerate(self.radiative.transitions)
                    if self.radiative.upperid[i] != id
                    and self.radiative.lowerid[i] != id
                ]
                self.radiative.upperid, self.radiative.lowerid = (
                    [
                        uid
                        for uid, lid in zip(
                            self.radiative.upperid, self.radiative.lowerid
                        )
                        if uid != id and lid != id
                    ],
                    [
                        lid
                        for uid, lid in zip(
                            self.radiative.upperid, self.radiative.lowerid
                        )
                        if uid != id and lid != id
                    ],
                )

                qns_str = ", ".join(
                    ["{}={}".format(key, val) for key, val in s.QN.__dict__.items()]
                )
                print("Warning: No matching energy levels found for " + qns_str)

        # update the data
        self.radiative.stateid = rad_stateid_new
        self.radiative.states = rad_states_new
        self.radiative.nstates = len(self.radiative.states)
        self.radiative.ntrans = len(self.radiative.transitions)
        self.stateid = self.radiative.stateid
        self.states = self.radiative.states
        self.nstates = self.radiative.nstates

        assert (
            self.collisional.nstates == self.nstates
        ), "Number of states do not match between data"
        assert set(self.stateid) == set(
            self.collisional.stateid
        ), "ID of states do not match between data"

    def write_LAMDA_file(self, filename, QNs=[]):
        # headers
        to_write = []
        to_write.append("!MOLECULE")
        to_write.append(self.formula)
        to_write.append("!MOLECULAR WEIGHT")
        to_write.append(
            str(self.molecular_weight) if self.molecular_weight is not None else ""
        )
        to_write.append("!NUMBER OF ENERGY LEVELS")
        to_write.append(str(self.nstates))
        line = "!LEVEL + ENERGY(CM-1) + WEIGHT + QUANTUM NOS."
        if len(QNs):
            line += "(" + ", ".join(QNs) + ")"
        to_write.append(line)

        # states
        for i, s in enumerate(self.states):
            qns_str = "_".join([str(s.QN.__dict__[key]) for key in QNs])
            line = "\t".join(
                [
                    str(i + 1).rjust(3),
                    "{:.6f}".format(float(s.energy)).rjust(11),
                    "{:.1f}".format(float(s.weight)).rjust(5),
                    qns_str,
                ]
            )
            to_write.append(line)

        # radiative transitions
        to_write.append("!NUMBER OF RADIATIVE TRANSITIONS")
        to_write.append(str(self.radiative.ntrans))
        to_write.append("!TRANS UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_up(K)")
        # sort
        upperidx = [self.stateid.index(id) for id in self.radiative.upperid]
        loweridx = [self.stateid.index(id) for id in self.radiative.lowerid]
        zipped = zip(
            upperidx,
            loweridx,
            self.radiative.transitions,
        )
        upperidx, loweridx, transitions = zip(*sorted(zipped))
        for i, trans in enumerate(transitions):
            uidx = upperidx[i] + 1
            lidx = loweridx[i] + 1
            line = "\t".join(
                [
                    str(i + 1).rjust(3),
                    str(uidx).rjust(3),
                    str(lidx).rjust(3),
                    "{:.4e}".format(trans.EinsteinA).rjust(12),
                    "{:.6f}".format(trans.restfreq * 1e-3).rjust(12),
                    "{:.3f}".format(trans.Eup).rjust(9),
                ]
            )
            to_write.append(line)

        # collisional transitions
        to_write.append("!NUMBER OF COLL PARTNERS")
        to_write.append(str(self.collisional.ncolpartners))

        for i, colpart in enumerate(list(set(self.collisional.colpartners))):
            # select only this collision partner
            indices = [
                i for i, cp in enumerate(self.collisional.colpartners) if cp == colpart
            ]
            transitions = [
                trans
                for i, trans in enumerate(self.collisional.transitions)
                if i in indices
            ]
            upperidx = [
                self.stateid.index(uid)
                for i, uid in enumerate(self.collisional.upperid)
                if i in indices
            ]
            loweridx = [
                self.stateid.index(lid)
                for i, lid in enumerate(self.collisional.lowerid)
                if i in indices
            ]
            # sort
            upperidx, loweridx, transitions = zip(
                *sorted(zip(upperidx, loweridx, transitions))
            )

            # comment and temperature
            comment = list(set([trans.comment for trans in transitions]))
            temperature = list(
                map(list, set(map(tuple, [trans.temperature for trans in transitions])))
            )
            assert (
                len(comment) == len(temperature) == 1
            ), "Comments or temperatures are not identical among the same collisional data"
            comment = comment[0]
            temperature = temperature[0]

            to_write.append("!COLLISIONS BETWEEN")
            to_write.append(comment)
            to_write.append("!NUMBER OF COLL TRANS")
            to_write.append(str(len(transitions)))
            to_write.append("!NUMBER OF COLL TEMPS")
            to_write.append(str(len(temperature)))
            to_write.append("!COLL TEMPS")
            to_write.append(
                "\t".join(["{:.1f}".format(t).rjust(6) for t in temperature])
            )
            to_write.append("!TRANS + UP + LOW + COLLRATES(cm^3 s^-1)")

            # coefficient
            lines = []
            for j, trans in enumerate(transitions):
                uidx = upperidx[j] + 1
                lidx = loweridx[j] + 1
                line = "\t".join(
                    [
                        str(j + 1).rjust(3),
                        str(uidx).rjust(3),
                        str(lidx).rjust(3),
                        "\t".join(
                            [
                                "{:.5e}".format(coeff).rjust(12)
                                for coeff in trans.coefficient
                            ]
                        ),
                    ]
                )
                lines.append(line)

            to_write.extend(lines)

        with open(filename, "w") as f:
            f.write("\n".join(to_write))

    def read_LAMDA_file(self, filename, QNs=[]):
        if not os.path.exists(filename):
            raise FileNotFoundError

        print("Reading the LAMDA file " + filename)
        f = open(filename, "r")
        # line 1-2; molecular formula
        f.readline()
        self.formula = f.readline()

        # line 3-4; MOLECULAR WEIGHT
        f.readline()
        self.molecular_weight = f.readline()

        # line 5-6; number of states (n)
        f.readline()
        self.nstates = int(f.readline())

        # line 7-7+n; states
        f.readline()  # header

        self.states = []
        self.stateid = []
        for _ in range(self.nstates):
            idx, energy, weight, qn_str = f.readline().strip().split()
            qn_list = qn_str.split("_")
            if len(QNs) != len(qn_list):
                raise ValueError(
                    "Mismatch between the numbers of provided QNs and QNs listed in the file"
                )
            qn_dict = {key: val for key, val in zip(QNs, qn_list)}
            QN = QuantumNumber(**qn_dict)
            state = State(energy=float(energy), weight=float(weight), QN=QN)
            self.states.append(state)
            self.stateid.append(idx)

        # line 7+n-7+n+1; number of radiative transitions (m)
        f.readline()
        ntrans = int(f.readline())

        # line 7+n+2-7+n+2+m; raditative transitions
        f.readline()
        lowerid = []
        upperid = []
        transitions = []
        for _ in range(ntrans):
            idx, uid, lid, EinsteinA, nu0, Eup = f.readline().strip().split()
            upperstate = self.states[self.stateid.index(uid)]
            lowerstate = self.states[self.stateid.index(lid)]
            trans = RadiativeTransition(
                species=self.formula,
                restfreq=float(nu0) * 1e9,  # in Hz
                restfreq_err=None,
                gup=upperstate.weight,
                glow=lowerstate.weight,
                EinsteinA=float(EinsteinA),
                Eup=float(Eup),
                QNup=upperstate.QN,
                QNlow=lowerstate.QN,
            )

            lowerid.append(lid)
            upperid.append(uid)
            transitions.append(trans)

        self.radiative = Transitions(
            transitions=transitions,
            states=self.states,
            lowerid=lowerid,
            upperid=upperid,
        )

        # line 7+n+2+m-7+n+2+m+1; number of collisional partners
        f.readline()
        ncolpartners = int(f.readline())

        colpartners = []
        lowerid = []
        upperid = []
        transitions = []
        for colpartID in range(ncolpartners):
            # description about the collisional partner
            f.readline()
            desc = f.readline()

            # number of transitions
            f.readline()
            ntrans = int(f.readline())

            # number of temperatures
            f.readline()
            ntemp = int(f.readline())

            # temperatures
            f.readline()
            temperatures = list(map(float, f.readline().strip().split()))

            # collisional rate
            f.readline()
            for _ in range(ntrans):
                data = f.readline().strip().split()
                idx, uid, lid = data[0:3]
                rates = list(map(float, data[3:]))

                assert ntemp == len(
                    rates
                ), "Mismatch between the numbers of temperature and collisional rate"
                upperstate = self.states[self.stateid.index(uid)]
                lowerstate = self.states[self.stateid.index(lid)]

                trans = CollisionalTransition(
                    species=self.formula,
                    temperature=temperatures,
                    coefficient=rates,
                    QNup=upperstate.QN,
                    QNlow=lowerstate.QN,
                    comment=desc,
                )

                colpartners.append(colpartID)
                lowerid.append(lid)
                upperid.append(uid)
                transitions.append(trans)

        print("Done.")

    # def select_data(self, restfreq=None, Eup=None, QN=None):
    #     # selection based on rest frequencies
    #     if restfreq is not None:
    #         restfreqs = [trans.restfreq for trans in self.radiative.transitions]
    #         numin, numax = restfreq
    #         indices = [i for i, nu in enumerate(restfreqs) if (nu >= numin) and (nu <= numax)]
