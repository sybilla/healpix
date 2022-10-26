using System;
using System.Collections.Generic;

namespace Sybilla.Healpix
{
    public static class HealpixExtensions
    {
        /// <summary>
        /// Query for pixels within disk (inclusive) given center provided as right ascension and declination angles and radius
        /// </summary>
        /// <param name="ra">right ascension in degrees</param>
        /// <param name="dec">declination in degrees</param>
        /// <param name="radius">radius in degrees</param>
        /// <returns></returns>
        public static List<int> QueryNested(int nside, double ra, double dec, double radius) 
        {
            var vec = Healpix.ang2vec(dec * Math.PI / 180 - Math.PI / 2, ra * Math.PI / 180);
            var list = new List<int>();
            var r = radius / 180 * Math.PI;
            Healpix.query_disc_inclusive_nest(nside, vec, r, list.Add);
            return list;
        }
    }
}
