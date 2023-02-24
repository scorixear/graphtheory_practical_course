import json

def write_json(obj, output):
    with open(output, "w", encoding="UTF-8") as writer:
        json.dump(obj, writer)
def read_json(input):
    with open(input, "r", encoding="UTF-8") as reader:
        return json.load(reader)
