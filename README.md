## Atmospheric electric satellite. See [<b>atmosat.ipynb</b>](https://github.com/gusgordon/atmosat/blob/master/atmosat.ipynb).

The goal of this project is to determine what the smallest solar plane that can indefinitely sustain flight is.

The approach is to first make a model of a solar plane that can be fed different parameters such as wingspan, battery mass, cruising altitude, number of props, payload mass, and so on. The model outputs the performance of the solar aircraft.

Different plane configutations are then fed into the model under an optimization procedure (differential evolution seems to work best). The goal of the optimization procedure is to find plane with the lowest wingspan that is able to indefinitely remain in the air.

### Results - see [notebook](https://github.com/gusgordon/atmosat/blob/master/atmosat.ipynb) for details. All units are in metric (meters, kg, kW, and so on).

![example_result](example_result.png)
