import yaml

def generate_yaml(yaml_path, data):
    with open(yaml_path, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)
