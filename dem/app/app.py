# Custom imports
from dem.src.core.factory import *
from dem.app.gravity      import *
from dem.app.restitution  import *
from dem.app.drop         import *
from dem.app.carreau      import *

# Declare factory
app_factory = factory()

# Register apps
app_factory.register("gravity",     gravity)
app_factory.register("restitution", restitution)
app_factory.register("drop",        drop)
app_factory.register("carreau",     carreau)
