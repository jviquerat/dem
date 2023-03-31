# Custom imports
from dem.src.core.factory     import *
from dem.src.domain.rectangle import *

# Declare factory
domain_factory = factory()

# Register domain
domain_factory.register("rectangle", rectangle)
