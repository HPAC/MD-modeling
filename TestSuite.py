from  General.Simulation import Simulation

ts = 1000
acc = 1e-4
mix = "geom"

TestSuite = {
    # Test SmallCube
    "SmallCube" :
        Simulation(
            "SmallCube",
            npart = 1000,
            domain = (11.01, 11.01, 11.01),
            box = (11.01, 11.01, 11.01),
            mixing = mix,
            restart = "config/e0_1000.rs",
            timesteps = 1000,
            nnodes = 1,
            nprocesses = 1,
            accuracy = acc
        ),
    # Test MediumCube
    "MediumCube" :
        Simulation(
            "MediumCube",
            npart = 64000,
            domain = (44.04, 44.04, 44.04),
            box = (44.04, 44.04, 44.04),
            mixing = mix, 
            restart = "config/e0_64000.rs",
            timesteps = ts,
            nnodes = 8,
            nprocesses = 64,
            accuracy = acc
        ),
    # Test LargeCube
    "LargeCube" :
        Simulation(
            "LargeCube",
            npart = 512000,
            domain = (88.08, 88.08, 88.08),
            box = (88.08, 88.08, 88.08),
            mixing = mix,
            restart = "config/e0_512000.rs",
            timesteps = ts,
            nnodes = 12,
            nprocesses = 96,
            accuracy = acc
        ),
    # Test LargeCube Sat
    "LargeCubeSat" :
        Simulation(
            "LargeCubeSat",
            npart = 512000,
            domain = (88.08, 88.08, 88.08),
            box = (88.08, 88.08, 88.08),
            mixing = mix,
            restart = "config/e0_512000.rs",
            timesteps = ts,
            nnodes = 3,
            nprocesses = 24,
            accuracy = acc
        ),
    # Test SmallSlab
    "SmallBulk" :
        Simulation(
            "SmallBulk",
            npart = 6000,
            domain = (11.01, 11.01, 66.06),
            box = (11.01, 11.01, 66.06),
            mixing = mix,
            restart = "config/e0_6000.rs",
            timesteps = ts,
            nnodes = 1,
            nprocesses = 8,
            accuracy = acc
        ),
    # Test MediumBulk
    "MediumBulk" :
        Simulation(
            "MediumBulk",
            npart = 64000,
            domain = (22.02, 22.02, 176.16),
            box = (22.02, 22.02, 176.16),
            mixing = mix,
            restart = "config/e0_64000.rs",
            timesteps = ts,
            nnodes = 1,
            nprocesses = 48,
            accuracy = acc
        ),
    # Test LargeSlab
    "LargeBulk" :
        Simulation(
            "LargeBulk",
            npart = 125000,
            domain = (44.04, 44.04, 176.16),
            box = (44.04, 44.04, 176.16),
            mixing = mix,
            restart = "config/e0_125000.rs",
            timesteps = ts,
            nnodes = 12,
            nprocesses = 96,
            accuracy = acc
        ),
    # Test SmallInterface
    "SmallInterface" :
        Simulation(
            "SmallInterface",
            npart = 4000,
            domain = (11.01, 11.01, 176.16),
            box = (11.01, 11.01, 44.04),
            mixing = mix,
            restart = "config/e0_4000.rs",
            timesteps = ts,
            nnodes = 1,
            nprocesses = 8,
            accuracy = 5e-4 #acc
        ),
    # Test MediumInterface
    "MediumInterface" :
        Simulation(
            "MediumInterface",
            npart = 32000,
            domain = (22.02, 22.02, 176.16),
            box = (22.02, 22.02, 88.08),
            mixing = mix,
            restart = "config/e0_32000.rs",
            timesteps = ts,
            nnodes = 8,
            nprocesses = 32,
            accuracy = acc
        ),
    # Case Study 1
    "CSBulk" :
        Simulation(
            "CSBulk",
            npart = 256000,
            domain = (55.00, 55.00, 110.00),
            box = (55.00, 55.00, 110.00),
            mixing = mix,
            restart = "config/e0_256000.rs",
            timesteps = ts,
            nnodes = 12,
            nprocesses = 96,
            accuracy = None,
            accuracy_real = 0.001,
            accuracy_kspace = 0.1
        )
}
