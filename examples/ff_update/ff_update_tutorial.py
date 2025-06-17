from pylift import utilities, reader

# reading in the gaff2 forcefield file packaged with PyLIFT
gaff2 = utilities.read_json("gaff2.json")

# reading in the gaff2 forcefield available directly in Amber and saving it in JSON format
gaff2 = reader.read_gaff2("gaff2.dat", default_loc=True, out_json="gaff2_updated.json")
