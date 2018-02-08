import os
if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
    import django
    django.setup()
    from viewer.loaders import load_from_dir
    load_from_dir("DCP2B","/code/media/DCP2B")
    load_from_dir("MURD", "/code/media/MURD")