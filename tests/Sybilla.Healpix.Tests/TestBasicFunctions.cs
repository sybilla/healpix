using System;
using System.Collections.Generic;
using Xunit;

namespace Sybilla.Healpix.Tests
{
    public class TestBasicFunctions
    {
        [Fact]
        public void test_nside_to_pixel_area()
        {
            var resolution = Healpix.nside2pixarea(256);
            Assert.Equal(1.5978966540475428e-05, resolution);
        }

        [Fact]
        public void test_nside_to_pixel_resolution()
        {
            var resolution = Healpix.nside2resol(256) * 180 * 60 / Math.PI;
            Assert.Equal(13.741945647269624, resolution,14);
        }

        [Fact]
        public void test_nside_to_npix()
        {
            var npix = Healpix.nside2npix(4);
            Assert.Equal(192, npix);
        }


        [Fact]
        public void test_cone_search_lonlat()
        {
            var pixels = HealpixExtensions.QueryNested(128, 45, 45, 1);
            var expectedPixels = new List<int>
            {
                12344,
                12339,
                12340,
                12333,
                12338,
                12337,
                12318,
                12332,
                12327,
                12336,
                12315,
                12316,
                12329,
                12326,
                12325,
                12314,
                12313,
                12310,
                12323,
                12324,
                12303,
                12312,
                12307,
                12321,
                12302,
                12301,
                12306,
                12299,
                12300,
                12295
            };
            expectedPixels.Sort();
            pixels.Sort();
            Assert.Equal(expectedPixels, pixels);
        }
    }
}
