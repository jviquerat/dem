# Custom imports
from dem.src.core.factory   import *
from dem.src.material.steel import *
from dem.src.material.glass import *

# Declare factory
material_factory = factory()

# Register materials
material_factory.register("steel", steel)
material_factory.register("glass", glass)
