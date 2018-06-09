import os


def get_target_list():
    if "TARGET_LIST" in os.environ:
        return os.environ["TARGET_LIST"].split(",")
    else:
        return ["CAMK1DA", "MURD", "HAO1A", "smTGR", "PTP1B"]
        # ,"DCLRE1AA","DCP2B","FALZA","MURD","NUDT21A",
        #    "NUDT22A","NUDT7A","PARP14A","PHIPA","SETDB1","SHMT2A","HAO1A","PTP1B",]


if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
    import django

    django.setup()
    from loader.loaders import process_target

    targets_to_load = get_target_list()
    prefix = "/code/media/"
    for target_name in targets_to_load:
        process_target(prefix, target_name)
