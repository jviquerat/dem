# Custom imports
from dem.src.core.factory import *
from dem.src.app.gravity  import *

# Declare factory
app_factory = factory()

# Register apps
app_factory.register("gravity", gravity)
