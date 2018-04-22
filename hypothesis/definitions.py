class IntTypes(object):

    def __init__(self):
        # The interaction class
        DEFAULT = 'DE'
        PROASIS = 'PR'
        self.INT_VER_CHOICES = (
                (DEFAULT, 'Default'),
                (PROASIS, 'Proasis'),
            )
        UNKNOWN = 'UK'
        HBOND = 'HB'
        VANDERWAALS = 'VD'
        CHARGE = 'CH'
        WEAKHBOND = 'WB'
        CATDIP = 'CD'
        CATPI = 'CP'
        DIPOLAR = 'DI'
        HALOGEN = 'HA'
        HDONPI = 'HP'
        PIPI = 'PP'
        POORANG = 'PA'
        UNFAV = 'UF'
        self.INT_TYPE_CHOICES = (
            (UNKNOWN, 'Unkwnon'),
            (HBOND, 'H-bond'),
            (CATPI, 'Cation Pi'),
            (HALOGEN, 'Halogen')
            (VANDERWAALS, 'Van der Waals'),
            (CHARGE, 'Charge'),
            (WEAKHBOND, 'Weak H-bond'),
            (CATDIP, 'Cation Dipole'),
            (HDONPI, 'H-bond donor Pi'),
            (DIPOLAR, 'Dipolar'),
            (PIPI, 'Pi Pi'),
            (POORANG, 'Poor angle'),
            (UNFAV, 'Unfavourable'),
        )
        self.DEFAULT_INT_VER = DEFAULT
        self.DEFAULT_INT_TYPE = UNKNOWN

        self.conv_dict = {PROASIS: {'cat_dip': CATDIP,
                                    'cat_pi': CATPI,
                                    'dipolar': DIPOLAR,
                                    'hbond': HBOND,
                                    'hdon_pi': HDONPI,
                                    'halogen': HALOGEN,
                                    'ionic': CHARGE,
                                    'pi_pi': PIPI,
                                    'poor_ang': POORANG,
                                    'unfav': UNFAV,
                                    'unk': UNKNOWN,
                                    'vdW': VANDERWAALS
                                    },
                          DEFAULT:{'cat_dip': CATDIP,
                                   'cat_pi': CATPI,
                                   'dipolar': DIPOLAR,
                                   'hbond': HBOND,
                                   'hdon_pi': HDONPI,
                                   'halogen': HALOGEN,
                                   'ionic': CHARGE,
                                   'pi_pi': PIPI,
                                   'poor_ang': POORANG,
                                   'unfav': UNFAV,
                                   'unk': UNKNOWN,
                                   'vdW': VANDERWAALS
                                   }
                          }

    def define_int_types(self):
        return self.INT_VER_CHOICES, self.DEFAULT_INT_VER, self.INT_TYPE_CHOICES, self.DEFAULT_INT_TYPE

    def get_int_conv(self,int_ver,int_name):
        if int_ver not in self.conv_dict:
            return None
        else:
            return self.conv_dict[int_ver][int_name]


class VectTypes(object):

    def __init__(self):
        ADDITION = 'AD'
        DELETION = 'DE'
        LINKER = 'LI'
        RING = 'RI'
        self.VECTOR_TYPES = (
            (ADDITION, 'Addition'),
            (DELETION, 'Deletion'),
            (LINKER, 'Linker'),
            (RING, 'Ring'),
        )
        self.vect_dict = {'additions': ADDITION, 'deletions': DELETION, 'linkers': LINKER, 'ring': RING}

    def get_vect_types(self):
        return self.VECTOR_TYPES

    def translate_vect_types(self,vect_input):
        return self.vect_dict[vect_input]
