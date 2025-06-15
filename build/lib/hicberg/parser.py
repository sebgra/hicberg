import configparser
from pathlib import Path

def parse_ini_file(file_path: Path) -> dict:
    """
    Parses an INI file and returns its content as a dictionary.
    """
    config = configparser.ConfigParser()
    try:
        config.read(file_path)
    except configparser.Error as e:
        print(f"Error parsing config file {file_path}: {e}")
        return {}

    parsed_data = {}
    for section in config.sections():
        parsed_data[section] = {}
        for key, value in config.items(section):
            # Attempt to convert types where appropriate
            if value.lower() in ('true', 'false'):
                parsed_data[section][key] = config.getboolean(section, key)
            elif value.isdigit():
                parsed_data[section][key] = config.getint(section, key)
            elif len(value.split(',')) != 1:
                parsed_data[section][key] = [e.strip() for e in value.split(',') if e.strip()]
            elif value == 'None' or value =='':
                parsed_data[section][key] = None
            else:
                try:
                    parsed_data[section][key] = config.getfloat(section, key)
                except ValueError:
                    parsed_data[section][key] = value # Keep as string if not other type

    return parsed_data