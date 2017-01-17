import unittest
import lsdmapartist

class TestLSDMapArtist(unittest.TestCase):

  def test_BaseRaster(self):
    self.assertIsInstance(lsdmapartist.BaseRaster)


