def calc1():
    from vtsim import vtsim as vt
    input = vt.read_json('./tutorial/tutorial_01.json')
    vt.run_calc(input)