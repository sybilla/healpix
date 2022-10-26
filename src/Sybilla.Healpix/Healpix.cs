using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Sybilla.Healpix
{
    // based on https://github.com/michitaro/healpix/blob/master/src/index.ts
    public static class Healpix
    {
        /**
         * # API Reference
         * 
         * This package based on this paper: [Gorski (2005)](http://iopscience.iop.org/article/10.1086/427976/pdf).
         *
         * The key things to understand the implementation are:
         * - Spherical coordinates in different representations such as `(alpha, delta)`
         *   or `(theta, phi)` or `(X, Y, z)` are always normalised to `(z, a)`.
         * - The HEALPix spherical projection is used to map to `(t, u)` (see `za2tu` and `tu2za`).
         *   See Section 4.4 and Figure 5 in the paper, where `(t, u)` is called `(x_s, y_s)`.
         *   
         * - A simple affine transformation is used to map to `(f, x, y)` (see `tu2fxy` and `fxy2tu`),
         *   where `f = {0 .. 11}` is the base pixel index and `(x, y)` is the position
         *   within the base pixel in the (north-east, north-west) direction
         *   and `(0, 0)` in the south corner.
         * - From `(f, x, y)`, the HEALPix pixel index in the "nested" scheme
         *   is related via `fxy2nest` and `nest2fxy`, and in the "ring" scheme
         *   via `fxy2ring` and `ring2fxy` in a relatively simple equations.
         * 
         * To summarise: there are two geometrical transformations:
         * `(z, a)` <-> `(t, u)` is the HEALPix spherical projection,
         * and `(t, u)` <-> `(f, x, y)` is a 45 deg rotation and scaling for each
         * of the 12 base pixels, so that HEALPix pixels in `(x, y)` are unit squares,
         * and pixel index compuatations are relatively straightforward,
         * both in the "nested" and "ring" pixelisation scheme.
         * 
         * ## Notations
         * 
         * <pre>
         * theta :  colatitude (pi/2 - delta)                [0 , pi]
         * phi   :  longitude (alpha)                        [0, 2 pi)
         * t     :  coord. of x-axis in spherical projection [0, 2 pi)
         * u     :  coord. of y-axis in spherical projection [-pi/2, pi/2]
         * z     :  cos(theta)                               [-1, 1]
         * X     :  sin(theta) * cos(phi)                    [-1, 1]
         * Y     :  sin(theta) * sin(phi)                    [-1, 1]
         * a     :  phi                                      [0, 2 pi)
         * f     :  base pixel index                         {0 .. 11}
         * x     :  north-east index in base pixel           [0, nside)
         * y     :  north-west index in base pixel           [0, nside)
         * p     :  north-east axis in base pixel            [0, 1)
         * q     :  north-west axis in base pixel            [0, 1)
         * j     :  pixel-in-ring index                      polar cap: {1 .. 4 i}
         *                                                   equatorial belt: {1 .. 4 nside}
         * i     :  ring index                               {1 .. 4 nside - 1}
         * </pre>
         */
        public static int order2nside(int order)
        {
            return 1 << order;
        }
        public static int nside2order(int nside)
        {
            return ilog2(nside);
        }
        public static int nside2npix(int nside)
        {
            return 12 * nside * nside;
        }
        public static int vec2pix_nest(int nside, double[] v)
        {
            var (z, a) = vec2za(v[0], v[1], v[2]);
            return za2pix_nest(nside, z, a);
        }
        public static int vec2pix_ring(int nside, double[] v)
        {
            var (z, a) = vec2za(v[0], v[1], v[2]);
            return nest2ring(nside, za2pix_nest(nside, z, a));
        }
        public static int ang2pix_nest(int nside, double theta, double phi)
        {
            var z = Math.Cos(theta);
            return za2pix_nest(nside, z, phi);
        }
        public static int ang2pix_ring(int nside, double theta, double phi)
        {
            var z = Math.Cos(theta);
            return nest2ring(nside, za2pix_nest(nside, z, phi));
        }
        public static int nest2ring(int nside, int ipix)
        {
            var (f, x, y) = nest2fxy(nside, ipix);
            return fxy2ring(nside, f, x, y);
        }
        public static int ring2nest(int nside, int ipix)
        {
            if (nside == 1)
            {
                return ipix;
            }
            var (f, x, y) = ring2fxy(nside, ipix);
            return fxy2nest(nside, f, x, y);
        }
        public static (int f, int x, int y) ring2fxy(int nside, int ipix)
        {
            var polar_lim = 2 * nside * (nside - 1);
            if (ipix < polar_lim)
            { // north polar cap
                var i = (int)Math.Floor((Math.Sqrt(1 + 2 * ipix) + 1) / 2);
                var j = ipix - 2 * i * (i - 1);
                var f = (int)Math.Floor((double)j / i);
                var k = j % i;
                var x = nside - i + k;
                var y = nside - 1 - k;
                return (f, x, y);
            }

            if (ipix < polar_lim + 8 * nside * nside)
            { // equatorial belt
                var k = ipix - polar_lim;
                var ring = 4 * nside;
                var i = nside - (int)Math.Floor((double)k / ring);
                var s = i % 2 == 0 ? 1 : 0;
                var j = 2 * (k % ring) + s;
                var jj = j - 4 * nside;
                var ii = i + 5 * nside - 1;
                var pp = (ii + jj) / 2;
                var qq = (ii - jj) / 2;
                var PP = (int)Math.Floor((double)pp / nside);
                var QQ = (int)Math.Floor((double)qq / nside);
                var V = 5 - (PP + QQ);
                var H = PP - QQ + 4;
                var f = 4 * V + (H >> 1) % 4;
                var x = pp % nside;
                var y = qq % nside;
                return (f, x, y);
            }
            else
            { // south polar cap
                var p = 12 * nside * nside - ipix - 1;
                var i = (int)Math.Floor((Math.Sqrt(1 + 2 * p) + 1) / 2);
                var j = p - 2 * i * (i - 1);
                var f = 11 - (int)Math.Floor((double)j / i);
                var k = j % i;
                var x = i - k - 1;
                var y = k;
                return (f, x, y);
            }
        }
        public static double[] pix2vec_nest(int nside, int ipix)
        {
            var (f, x, y) = nest2fxy(nside, ipix);
            var (t, u) = fxy2tu(nside, f, x, y);
            var (z, a) = tu2za(t, u);
            return za2vec(z, a);
        }
        public static (double theta, double phi) pix2ang_nest(int nside, int ipix)
        {
            var (f, x, y) = nest2fxy(nside, ipix);
            var (t, u) = fxy2tu(nside, f, x, y);
            var (z, a) = tu2za(t, u);
            return (theta: Math.Acos(z), phi: a);
        }
        public static double[] pix2vec_ring(int nside, int ipix)
        {
            return pix2vec_nest(nside, ring2nest(nside, ipix));
        }
        public static (double theta, double phi) pix2ang_ring(int nside, int ipix)
        {
            return pix2ang_nest(nside, ring2nest(nside, ipix));
        }
        // TODO: cleanup    
        public static void query_disc_inclusive_nest(int nside, double[] v, double radius, Action<int> cb)
        {
            if (radius > PI_2)
            {
                throw new Exception("query_disc: radius must < PI/2");
            }
            var pixrad = max_pixrad(nside);
            var d = PI_4 / nside;
            var (z0, a0) = vec2za(v[0], v[1], v[2]); // z0 = cos(theta)
            var sin_t = Math.Sqrt(1 - z0 * z0);
            var cos_r = Math.Cos(radius); // r := radius
            var sin_r = Math.Sin(radius);
            var z1 = z0 * cos_r + sin_t * sin_r; // cos(theta - r)
            var z2 = z0 * cos_r - sin_t * sin_r; // cos(theta + r)
            var u1 = za2tu(z1, 0).u;
            var u2 = za2tu(z2, 0).u;
            var cover_north_pole = sin_t * cos_r - z0 * sin_r < 0; // sin(theta - r) < 0
            var cover_south_pole = sin_t * cos_r + z0 * sin_r < 0; // sin(theta - r) < 0
            var i1 = (int)Math.Floor((PI_2 - u1) / d);
            var i2 = (int)Math.Floor((PI_2 - u2) / d + 1);
            if (cover_north_pole)
            {
                ++i1;
                for (var i = 1; i <= i1; ++i)
                    walk_ring(nside, i, cb);
                ++i1;
            }
            if (i1 == 0)
            {
                walk_ring(nside, 1, cb);
                i1 = 2;
            }
            if (cover_south_pole)
            {
                --i2;
                for (var i = i2; i <= 4 * nside - 1; ++i)
                    walk_ring(nside, i, cb);
                --i2;
            }
            if (i2 == 4 * nside)
            {
                walk_ring(nside, 4 * nside - 1, cb);
                i2 = 4 * nside - 2;
            }
            var theta = Math.Acos(z0);
            for (var i = i1; i <= i2; ++i)
                walk_ring_around(nside, i, a0, theta, radius + pixrad, ipix =>
                {
                    var temp = pix2vec_nest(nside, ipix);
                    var ang = angle(temp, v);
                    if (ang <= radius + pixrad)
                        cb(ipix);
                });
        }
        public static void query_disc_inclusive_ring(int nside, double[] v, double radius, Action<int> cb_ring)
        {
            query_disc_inclusive_nest(nside, v, radius, ipix =>
            {
                cb_ring(nest2ring(nside, ipix));
            });
        }
        public static double max_pixrad(int nside)
        {
            var unit = PI_4 / nside;
            var a = tu2vec(unit, nside * unit);
            var b = tu2vec(unit, (nside + 1) * unit);
            return angle(a
                ,
                b
             );
        }
        private static double angle(double[] a, double[] b)
        {
            return 2 * Math.Asin(Math.Sqrt(distance2(a, b)) / 2);
        }
        private static double[] tu2vec(double t, double u)
        {
            var (z, a) = tu2za(t, u);
            return za2vec(z, a);
        }
        private static double distance2(double[] a, double[] b)
        {
            var dx = a[0] - b[0];
            var dy = a[1] - b[1];
            var dz = a[2] - b[2];
            return dx * dx + dy * dy + dz * dz;
        }
        private static void walk_ring_around(int nside, double i, double a0, double theta, double r, Action<int> cb)
        {
            if (theta < r || theta + r > PI)
            {
                walk_ring(nside, i, cb);
                return;
            }

            var u = PI_4 * (2 - i / nside);
            var z = tu2za(PI_4, u).z;
            var st = Math.Sin(theta);
            var ct = Math.Cos(theta);
            var sr = Math.Sin(r);
            var cr = Math.Cos(r);
            var w = Math.Atan2(
                Math.Sqrt(-square(z - ct * cr) / (square(st) * sr * sr) + 1) * sr,
                (-z * ct + cr) / st
            );
            if (w >= PI)
            {
                walk_ring(nside, i, cb);
                return;
            }
            var t1 = center_t(nside, i, za2tu(z, wrap(a0 - w, PI2)).t);
            var t2 = center_t(nside, i, za2tu(z, wrap(a0 + w, PI2)).t);
            var begin = tu2fxy(nside, t1, u);
            var temp = tu2fxy(nside, t2, u);
            var end = right_next_pixel(nside, temp.f, temp.x, temp.y);
            for (var s = begin; !fxy_compare(s, end); s = right_next_pixel(nside, s.f, s.x, s.y))
            {
                var nest = fxy2nest(nside, s.f, s.x, s.y);
                cb(nest);
            }
        }
        private static double center_t(int nside, double i, double t)
        {
            var d = PI_4 / nside;
            t /= d;
            t = (((int)(t + i % 2) >> 1) << 1) + 1 - i % 2; // potential bug here
            t *= d;
            return t;
        }
        private static void walk_ring(int nside, double i, Action<int> cb)
        {
            var u = PI_4 * (2 - i / nside);
            var t = PI_4 * (1 + (1 - i % 2) / nside);
            var begin = tu2fxy(nside, t, u);
            var s = begin;
            do
            {
                cb(fxy2nest(nside, s.f, s.x, s.y));
                s = right_next_pixel(nside, s.f, s.x, s.y);
            } while (!fxy_compare(s, begin));
        }
        private static bool fxy_compare((int f, int x, int y) a, (int f, int x, int y) b)
        {
            return a.x == b.x && a.y == b.y && a.f == b.f;
        }
        private static (int f, int x, int y) right_next_pixel(int nside, int f, int x, int y)
        {
            ++x;
            if (x == nside)
            {
                switch ((int)Math.Floor((double)f / 4))
                {
                    case 0:
                        f = (f + 1) % 4;
                        x = y;
                        y = nside;
                        break;
                    case 1:
                        f = f - 4;
                        x = 0;
                        break;
                    case 2:
                        f = 4 + (f + 1) % 4;
                        x = 0;
                        break;
                }
            }
            --y;
            if (y == -1)
            {
                switch ((int)Math.Floor((double)f / 4))
                {
                    case 0:
                        f = 4 + (f + 1) % 4;
                        y = nside - 1;
                        break;
                    case 1:
                        f = f + 4;
                        y = nside - 1;
                        break;
                    case 2:
                        {
                            f = 8 + (f + 1) % 4;
                            y = x - 1;
                            x = 0;
                            break;
                        }
                }
            }
            return (f, x, y);
        }
        public static IReadOnlyList<double[]> corners_nest(int nside, int ipix)
        {
            var (f, x, y) = nest2fxy(nside, ipix);
            var (t, u) = fxy2tu(nside, f, x, y);
            var d = PI_4 / nside;
            var xyzs = new List<double[]>();
            var corners = new List<(double tt, double uu)>
            {
                (0, d),
                (-d, 0),
                (0, -d),
                (d, 0)
            };

            foreach (var (tt, uu) in corners)
            {
                var (z, a) = tu2za(t + tt, u + uu);
                xyzs.Add(za2vec(z, a));
            }
            return xyzs;
        }
        public static IReadOnlyList<double[]> corners_ring(int nside, int ipix)
        {
            return corners_nest(nside, ring2nest(nside, ipix));
        }
        // pixel area
        public static double nside2pixarea(int nside)
        {
            return PI / (3 * nside * nside);
        }
        // average pixel size
        public static double nside2resol(int nside)
        {
            return Math.Sqrt(PI / 3) / nside;
        }
        public static double[] pixcoord2vec_nest(int nside, int ipix, int ne, int nw)
        {
            var (f, x, y) = nest2fxy(nside, ipix);
            var (t, u) = fxy2tu(nside, f, x, y);
            var d = PI_4 / nside;
            var (z, a) = tu2za(t + d * (ne - nw), u + d * (ne + nw - 1));
            return za2vec(z, a);
        }
        public static double[] pixcoord2vec_ring(int nside, int ipix, int ne, int nw)
        {
            return pixcoord2vec_nest(nside, ring2nest(nside, ipix), ne, nw);
        }
        private static int za2pix_nest(int nside, double z, double a)
        {
            var (t, u) = za2tu(z, a);
            var (f, x, y) = tu2fxy(nside, t, u);
            return fxy2nest(nside, f, x, y);
        }
        public static (int f, int x, int y) tu2fxy(int nside, double t, double u)
        {
            var (f, p, q) = tu2fpq(t, u);
            var x = (int)clip(Math.Floor(nside * p), 0, nside - 1);
            var y = (int)clip(Math.Floor(nside * q), 0, nside - 1);
            return (f, x, y);
        }
        private static double wrap(double A, double B)
        {
            return A < 0 ? B - (-A % B) : A % B;
        }

        private static double PI2 = 2 * Math.PI;
        private static double PI = Math.PI;
        private static double PI_2 = Math.PI / 2;
        private static double PI_4 = Math.PI / 4;
        private static double PI_8 = Math.PI / 8;

        private static double sigma(double z)
        {
            if (z < 0)
                return -sigma(-z);
            else
                return 2 - Math.Sqrt(3 * (1 - z));
        }

        /**
       * HEALPix spherical projection.
       */
        public static (double t, double u) za2tu(double z, double a)
        {
            if (Math.Abs(z) <= 2.0 / 3.0)
            { // equatorial belt
                var t = a;
                var u = 3 * PI_8 * z;
                return (t, u);
            }
            else
            { // polar caps
                var p_t = a % (PI_2);
                var sigma_z = sigma(z);
                var t = a - (Math.Abs(sigma_z) - 1) * (p_t - PI_4);
                var u = PI_4 * sigma_z;
                return (t, u);
            }
        }

        /**
         * Inverse HEALPix spherical projection.
         */
        public static (double z, double a) tu2za(double t, double u)
        {
            var abs_u = Math.Abs(u);
            if (abs_u >= PI_2)
            { // error
                return (z: sign(u), a: 0);
            }
            if (abs_u <= Math.PI / 4)
            { // equatorial belt
                var z = 8 / (3 * PI) * u;
                var a = t;
                return (z, a);
            }
            else
            { // polar caps
                var t_t = t % (Math.PI / 2);
                var a = t - (abs_u - PI_4) / (abs_u - PI_2) * (t_t - PI_4);
                var z = sign(u) * (1 - 1.0 / 3 * square(2 - 4 * abs_u / PI));
                return (z, a);
            }
        }

        // (x, y, z) -> (z = cos(theta), phi)
        public static (double z, double a) vec2za(double X, double Y, double z)
        {
            var r2 = X * X + Y * Y;
            if (r2 == 0)
                return (z: z < 0 ? -1 : 1, a: 0);
            else
            {
                var a = (Math.Atan2(Y, X) + PI2) % PI2;
                z /= Math.Sqrt(z * z + r2);
                return (z, a);
            }
        }

        // (z = cos(theta), phi) -> (x, y, z)
        public static double[] za2vec(double z, double a)
        {
            var sin_theta = Math.Sqrt(1 - z * z);
            var X = sin_theta * Math.Cos(a);
            var Y = sin_theta * Math.Sin(a);
            return new double[] { X, Y, z };
        }
        public static double[] ang2vec(double theta, double phi)
        {
            var z = Math.Cos(theta);
            return za2vec(z, phi);
        }

        public static (double theta, double phi) vec2ang(double[] v)
        {
            var (z, a) = vec2za(v[0], v[1], v[2]);
            return (theta: Math.Acos(z), phi: a);
        }

        // spherical projection -> f, p, q
        // f: base pixel index
        // p: coord in north east axis of base pixel
        // q: coord in north west axis of base pixel
        private static (int f, double p, double q) tu2fpq(double t, double u)
        {
            t /= PI_4;
            u /= PI_4;
            t = wrap(t, 8);
            t += -4;
            u += 5;
            var pp = clip((u + t) / 2, 0, 5);
            var PP = (int)Math.Floor(pp);
            var qq = clip((u - t) / 2, 3 - PP, 6 - PP);
            var QQ = (int)Math.Floor(qq);
            var V = 5 - (PP + QQ);
            if (V < 0)
            { // clip
                return (f: 0, p: 1, q: 1);
            }
            var H = PP - QQ + 4;
            var f = 4 * V + (H >> 1) % 4;
            var p = pp % 1;
            var q = qq % 1;
            return (f, p, q);
        }
        // f, p, q -> nest index
        public static int fxy2nest(int nside, int f, int x, int y)
        {
            return f * nside * nside + bit_combine(x, y);
        }
        // x = (...x2 x1 x0)_2 <- in binary
        // y = (...y2 y1 y0)_2
        // p = (...y2 x2 y1 x1 y0 x0)_2
        // returns p
        public static int bit_combine(int x, int y)
        {
            assert(x < (1 << 16));
            assert(y < (1 << 15));

            return (
                // (python)
                // n = 14
                // ' | '.join(['x & 1'] + [f'(x & 0x{2 ** (i+1):x} | y & 0x{2 ** i:x}) << {i + 1}' for i in range(n)] + [f'y & 0x{2**n:x} << {n+1}'])
                x & 1 | (x & 0x2 | y & 0x1) << 1 | (x & 0x4 | y & 0x2) << 2 |
                (x & 0x8 | y & 0x4) << 3 | (x & 0x10 | y & 0x8) << 4 | (x & 0x20 | y & 0x10) << 5 |
                (x & 0x40 | y & 0x20) << 6 | (x & 0x80 | y & 0x40) << 7 | (x & 0x100 | y & 0x80) << 8 |
                (x & 0x200 | y & 0x100) << 9 | (x & 0x400 | y & 0x200) << 10 | (x & 0x800 | y & 0x400) << 11 |
                (x & 0x1000 | y & 0x800) << 12 | (x & 0x2000 | y & 0x1000) << 13 | (x & 0x4000 | y & 0x2000) << 14 |
                (x & 0x8000 | y & 0x4000) << 15 | y & 0x8000 << 16
            );
        }
        // x = (...x2 x1 x0)_2 <- in binary
        // y = (...y2 y1 y0)_2
        // p = (...y2 x2 y1 x1 y0 x0)_2
        // returns x, y
        public static (int x, int y) bit_decombine(int p)
        {
            //assert(p <= 0x7fffffff);
            // (python)
            // ' | '.join(f'(p & 0x{2**(2*i):x}) >> {i}' for i in range(16))
            var x = (p & 0x1) >> 0 | (p & 0x4) >> 1 | (p & 0x10) >> 2 |
        (p & 0x40) >> 3 | (p & 0x100) >> 4 | (p & 0x400) >> 5 |
        (p & 0x1000) >> 6 | (p & 0x4000) >> 7 | (p & 0x10000) >> 8 |
        (p & 0x40000) >> 9 | (p & 0x100000) >> 10 | (p & 0x400000) >> 11 |
        (p & 0x1000000) >> 12 | (p & 0x4000000) >> 13 | (p & 0x10000000) >> 14 | (p & 0x40000000) >> 15;
            // (python)
            // ' | '.join(f'(p & 0x{2**(2*i + 1):x}) >> {i+1}' for i in range(15))
            var y = (p & 0x2) >> 1 | (p & 0x8) >> 2 | (p & 0x20) >> 3 |
        (p & 0x80) >> 4 | (p & 0x200) >> 5 | (p & 0x800) >> 6 |
        (p & 0x2000) >> 7 | (p & 0x8000) >> 8 | (p & 0x20000) >> 9 |
        (p & 0x80000) >> 10 | (p & 0x200000) >> 11 | (p & 0x800000) >> 12 |
        (p & 0x2000000) >> 13 | (p & 0x8000000) >> 14 | (p & 0x20000000) >> 15;
            return (x, y);
        }
        // f: base pixel index
        // x: north east index in base pixel
        // y: north west index in base pixel
        private static (int f, int x, int y) nest2fxy(int nside, int ipix)
        {
            var nside2 = nside * nside;
            var f = (int)Math.Floor((double)ipix / nside2); // base pixel index
            var k = ipix % nside2;             // nested pixel index in base pixel
            var (x, y) = bit_decombine(k);
            return (f, x, y);
        }
        private static int fxy2ring(int nside, int f, int x, int y)
        {
            var f_row = (int)Math.Floor(f / 4.0); // {0 .. 2}
            var f1 = f_row + 2;            // {2 .. 4}
            var v = x + y;
            var i = f1 * nside - v - 1;

            if (i < nside)
            { // north polar cap
                var f_col = f % 4;
                var ipix = 2 * i * (i - 1) + (i * f_col) + nside - y - 1;
                return ipix;
            }
            if (i < 3 * nside)
            { // equatorial belt
                var h = x - y;
                var f2 = 2 * (f % 4) - (f_row % 2) + 1;  // {0 .. 7}
                var k = (f2 * nside + h + (8 * nside)) % (8 * nside);
                var offset = 2 * nside * (nside - 1);
                var ipix = offset + (i - nside) * 4 * nside + (k >> 1);
                return ipix;
            }
            else
            { // south polar cap
                var i_i = 4 * nside - i;
                var i_f_col = 3 - (f % 4);
                var j = 4 * i_i - (i_i * i_f_col) - y;
                var i_j = 4 * i_i - j + 1;
                var ipix = 12 * nside * nside - 2 * i_i * (i_i - 1) - i_j;
                return ipix;
            }
        }
        // f, x, y -> spherical projection
        public static (double t, double u) fxy2tu(int nside, int f, int x, int y)
        {
            var f_row = (int)Math.Floor((double)f / 4);
            var f1 = f_row + 2;
            var f2 = 2 * (f % 4) - (f_row % 2) + 1;
            var v = x + y;
            var h = x - y;
            var i = f1 * nside - v - 1;
            var k = (f2 * nside + h + (8 * nside));
            var t = (double)k / nside * PI_4;
            var u = PI_2 - (double)i / nside * PI_4;
            return (t, u);
        }
        public static long orderpix2uniq(int order, int ipix)
        {
            /**
             * Pack `(order, ipix)` into a `uniq` integer.
             * 
             * This HEALPix "unique identifier scheme" is starting to be used widely:
             * - see section 3.2 in http://healpix.sourceforge.net/pdf/intro.pdf
             * - see section 2.3.1 in http://ivoa.net/documents/MOC/
             */
            return 4 * ((1 << (2 * order)) - 1) + ipix;
        }
        public static (int order, int ipix) uniq2orderpix(int uniq)
        {
            /**
             * Unpack `uniq` integer into `(order, ipix)`.
             * 
             * Inverse of `orderpix2uniq`.
             */
            assert(uniq <= 0x7fffffff);
            var order = 0;
            var l = (uniq >> 2) + 1;
            while (l >= 4)
            {
                l >>= 2;
                ++order;
            }
            var ipix = uniq - (((1 << (2 * order)) - 1) << 2);
            return (order, ipix);
        }
        public static int ilog2(int x)
        {
            /**
             * log2 for integer numbers.
             * 
             * We're not calling Math.log2 because it's not supported on IE yet.
             */
            var o = -1;
            while (x > 0)
            {
                x >>= 1;
                ++o;
            }
            return o;
        }
        private static double sign(double A)
        {
            return A > 0 ? 1 : (A < 0 ? -1 : 0);
        }
        private static double square(double A)
        {
            return A * A;
        }
        private static double clip(double Z, double A, double B)
        {
            return Z < A ? A : (Z > B ? B : Z);
        }
        private static void assert(bool condition)
        {
            Debug.Assert(condition);
            if (!condition)
            {
                throw new Exception("Assertion failed");
            }
        }

    }

































}
