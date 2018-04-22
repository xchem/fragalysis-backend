import os

if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
    import django
    django.setup()
    from loader.loaders import load_from_dir,load_events_from_dir
    from loader.loaders import analyse_target
    targets_to_load = ["CAMK1DA","DCLRE1AA","DCP2B","FALZA","MURD","NUDT21A","NUDT22A","NUDT7A",
                       "PARP14A","PHIPA","SETDB1","SHMT2A","HAO1A","PTP1B",]
    for t in targets_to_load:
        load_from_dir(t,"/code/media/"+t)
        load_events_from_dir(t,"/code/media/"+t+"/events")

    for t in targets_to_load:
        analyse_target(t)
