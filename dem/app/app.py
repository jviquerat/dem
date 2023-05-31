# Custom imports
from dem.src.core.factory import *
from dem.app.gravity      import *
from dem.app.restitution  import *
from dem.app.dam_break    import *
from dem.app.carreau      import *

# Declare factory
app_factory = factory()

# Register apps
app_factory.register("gravity",     gravity)
app_factory.register("restitution", restitution)
app_factory.register("dam_break",   dam_break)
app_factory.register("carreau",     carreau)
