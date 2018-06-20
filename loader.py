import os


def get_target_list(base_path):
    target_path = os.path.join(base_path, "TARGET_LIST")
    if os.path.isfile(target_path):
        return [x.strip() for x in open(target_path).read().split(" ") if x.strip()]
    elif "TARGET_LIST" in os.environ:
        return os.environ["TARGET_LIST"].split(",")
    else:
        return ["MURD", "HAO1A", "smTGR", "PTP1B"]


if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
    import django

    django.setup()
    from loader.loaders import process_target

    prefix = "/code/media/"
    targets_to_load = get_target_list(prefix)
    for target_name in targets_to_load:
        process_target(prefix, target_name)
