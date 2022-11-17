# Custom imports
from dem.src.core.factory    import *
from dem.src.app.gravity     import *
from dem.src.app.restitution import *

# Declare factory
app_factory = factory()

# Register apps
app_factory.register("gravity",     gravity)
app_factory.register("restitution", restitution)
