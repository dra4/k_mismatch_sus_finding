import json


def load_json(fname):
    with open(fname) as f:
        data = f.read()
        return json.loads(data)


def load_fasta(fname):
    with open(fname) as f:
        f.readline()
        data = f.read().replace("\n", "")
        return data


def write_fasta(fname, sx):
    with open(fname, 'w') as f:
        f.write(">0\n")
        f.write(sx)
        f.write("\n")
