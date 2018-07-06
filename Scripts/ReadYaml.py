import yaml

with open("DataMCInfo.yml","r") as yamlfile:
	cfg = yaml.load(yamlfile)

for section in cfg:
    print(section)
    print(cfg[section])
