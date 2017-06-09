import mitsuba
from mitsuba.core import *
from mitsuba.render import SceneHandler
fileResolver = Thread.getThread().getFileResolver()
fileResolver.appendPath('~/CONF/QUIJIBO/renders/python')
paramMap = StringMap()
scene = SceneHandler.loadScene(fileResolver.resolve("/Users/tkim/CONF/QUIJIBO/renders/python/bunny_zoom_0_1.xml"), paramMap)

from mitsuba.core import *
from mitsuba.render import RenderQueue, RenderJob
import multiprocessing

scheduler = Scheduler.getInstance()
# Start up the scheduling system with one worker per local core
for i in range(0, 4):
    scheduler.registerWorker(LocalWorker(i, 'wrk%i' % i))
scheduler.start()
queue = RenderQueue()

import numpy as np
from scipy import *

target0 = Vector(0.681017,0.344339, 0.841651)
origin0 = Vector(1.06275, -0.147389, 1.62427)
up0 = Vector(-0.180728, -0.870101, -0.458544)

target1 = Vector(0.0329559,0.23162,0.731513) 
origin1 = Vector(0.975366,0.0865425,0.430146)
up1 = Vector(-0.0605996, -0.960187,0.272722)

sensor = scene.getSensor()

frames = 60
for x in range(0,frames):
    frac = float(x) / float(frames - 1)
    target = (1.0 - frac) * target0 + frac * target1
    origin = (1.0 - frac) * origin0 + frac * origin1
    up     = (1.0 - frac) * up0     + frac * up1
    xform = Transform.lookAt(Point(origin), Point(target), up)

    fracNext = float(x + 1) / float(frames - 1)
    targetNext = (1.0 - fracNext) * target0 + fracNext * target1
    originNext = (1.0 - fracNext) * origin0 + fracNext * origin1
    upNext     = (1.0 - fracNext) * up0     + fracNext * up1
    xformNext = Transform.lookAt(Point(originNext), Point(targetNext), upNext)

    animatedXform = AnimatedTransform()
    animatedXform.appendTransform(0, xform)
    animatedXform.appendTransform(1, xformNext)
    animatedXform.sortAndSimplify()

    sensor.setWorldTransform(animatedXform)

    # below here seems to be what needs to be in the inner loop
    filename = "result_motionblur_" + str(x) + ".exr"
    scene.setDestinationFile(filename)
    job = RenderJob('myRenderJob', scene, queue)
    job.start()
    queue.waitLeft(0)
    queue.join()
    print(Statistics.getInstance().getStats())
