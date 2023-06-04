# Custom imports
from dem.src.core.factory     import *
from dem.src.domain.rectangle import *
from dem.src.domain.circle    import *

# Declare factory
domain_factory = factory()

# Register domain
domain_factory.register("rectangle", rectangle)
domain_factory.register("circle",    circle)
