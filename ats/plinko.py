from yggdrasil import *

# Plinko: particles fall under constant gravity through a staggered maze of
# sphere/box colliders. Not analytically solvable (the bounce sequence is
# effectively chaotic), but fully deterministic, so this is a snapshot
# regression test -- it returns to the same state at the chosen stop time every
# run. (This is the collider-maze that used to live in kinetics_grav.py.)


def run():
    g = -9.8
    numNodes = 20
    constants = MKS()
    nodes = NodeList(numNodes)

    constantForce = ConstantForce2d(nodes, constants, Vector2d(0, g))
    kinetics = Kinetics2d(nodes, constants)
    packages = [constantForce, kinetics]

    spacing_x = 2
    spacing_y = (6 - 2) / 5
    rows = 8
    collider_radius = 0.3
    cbounds = []
    for i in range(rows):
        num = 10 if i % 2 == 0 else 9
        offset = 0 if i % 2 == 0 else spacing_x / 2
        for j in range(num):
            x = (j - num // 2) * spacing_x + offset
            y = 2 + i * spacing_y
            cbounds.append(SphereCollider2d(position=Vector2d(x, y), radius=collider_radius, elasticity=0.8))
    cbounds.append(BoxCollider2d(position1=Vector2d(-10, -3), position2=Vector2d(10.5, 1), elasticity=0.5))
    cbounds.append(BoxCollider2d(position1=Vector2d(-11, 0), position2=Vector2d(-9.5, 10), elasticity=1))
    cbounds.append(BoxCollider2d(position1=Vector2d(9.5, 0), position2=Vector2d(11, 10), elasticity=1))
    for bound in cbounds:
        kinetics.addBoundary(bound)

    integrator = RungeKutta2Integrator2d(packages=packages, dtmin=0.01, verbose=False)

    rad = nodes.getFieldDouble("radius")
    mass = nodes.getFieldDouble("mass")
    for i in range(numNodes):
        rad[i] = 0.2
        mass[i] = 0.2

    pos = nodes.getFieldVector2d("position")
    for i in range(numNodes):
        pos[i] = Vector2d(-9 + i * (18 / numNodes), 10)

    Controller(integrator=integrator, periodicWork=[], statStep=200).Step(150)

    values = [integrator.time]
    for i in range(numNodes):
        values.append(pos[i].y)
    return {"mode": "snapshot", "values": values}


if __name__ == "__main__":
    print(run())
