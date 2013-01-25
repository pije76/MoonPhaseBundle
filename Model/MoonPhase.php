<?php

namespace Swis\Bundle\MoonPhaseBundle\Model;

/**
 * Moon phase calculation class
 * Adapted for PHP from Moontool for Windows (http://www.fourmilab.ch/moontoolw/)
 * by Samir Shah (http://rayofsolaris.net)
 * Last modified August 2012
 * 
 * @link https://github.com/solarissmoke/php-moon-phase
 * @author Samir Shah (http://rayofsolaris.net)
 * @license MIT
 */
class MoonPhase
{
    /*
     * Astronomical constants
     */

    const EPOCH = 2444238.5;   // 1980 January 0.0

    /*
     * Constants defining the Sun's apparent orbit
     */
    const E_LONG_EPOCH = 278.833540;  // Ecliptic longitude of the Sun at epoch 1980.0
    const E_LONG_PERIGEE = 282.596403;  // Ecliptic longitude of the Sun at perigee
    const ECC = 0.016718;   // Eccentricity of Earth's orbit
    const SUNSMAX = 1.495985e8;  // Semi-major axis of Earth's orbit, km
    const SUNANGSIZ = 0.533128;  // Sun's angular size, degrees, at semi-major axis distance

    /*
     * Elements of the Moon's orbit, epoch 1980.0
     */
    const MMLONG = 64.975464;  // Moon's mean longitude at the epoch
    const MMLONGP = 349.383063;  // Mean longitude of the perigee at the epoch
    const MLNODE = 151.950429;  // Mean longitude of the node at the epoch
    const MINC = 5.145396;   // Inclination of the Moon's orbit
    const MECC = 0.054900;   // Eccentricity of the Moon's orbit
    const MANGSIZ = 0.5181;   // Moon's angular size at distance a from Earth
    const MSMAX = 384401;   // Semi-major axis of Moon's orbit in km
    const MPARALLAX = 0.9507;  // Parallax at distance a from Earth
    const SYNODIC_MONTH = 29.53058868; // Synodic month (new Moon to new Moon)
    const LUNATBASE = 2423436.0;  // Base date for E. W. Brown's numbered series of lunations (1923 January 16)    

    /**
     * Current time [s]
     *  
     * @var number
     */

    private $timestamp;

    /**
     * Moon phase (0 to 1)
     * 
     * @var number
     */
    private $phase;

    /**
     * Amount of illuminated surface (0 to 1)
     *  
     * @var number
     */
    private $illum;

    /**
     * Moon age [d]
     *  
     * @var number
     */
    private $age;

    /**
     * Distance to Earth [km]
     *  
     * @var number
     */
    private $dist;

    /**
     * Angular diameter [rad]
     *  
     * @var number
     */
    private $angdia;

    /**
     * Distance to Sun [km]
     *  
     * @var number
     */
    private $distanceToSun;

    /**
     * Angular diameter of the sun [rad]
     * 
     * @var number
     */
    private $sunAngularDiameter;

    /**
     * Lunar cycle timestamps for each quarter [s]
     *  
     * @var array
     */
    private $quarters = null;
    private $timezone = null;

    /**
     * Calculate all data for the current time (Or an optional timestamp)
     *  
     * @param number $date
     */
    function __construct(\DateTime $date = null)
    {
        if (\is_null($date)) {
            $date = new \DateTime();
        }



        /*
         * General calculations
         */

        $this->timezone = $date->getTimezone();
        $this->timestamp = $date->format('U');
        $dateWithinEpoch = $this->utcToJulian($this->timestamp) - self::EPOCH;

        /*
         * Calculation of the Sun's position
         */

        $meanAnomalyOfSun = $this->normalizeAngle((360 / 365.2422) * $dateWithinEpoch);
        $perigeeToEpoch = $this->normalizeAngle($meanAnomalyOfSun + self::E_LONG_EPOCH - self::E_LONG_PERIGEE);
        $keplerSolution = $this->solveKeplerEquation($perigeeToEpoch, self::ECC);
        $trueAnomalyHelper = \sqrt((1 + self::ECC) / (1 - self::ECC)) * \tan($keplerSolution / 2);
        $trueAnomaly = 2 * \rad2deg(\atan($trueAnomalyHelper));
        $eclipticLong = $this->normalizeAngle($trueAnomaly + self::E_LONG_PERIGEE);

        $orbitalDistanceFactor = ((1 + self::ECC * cos(deg2rad($trueAnomaly))) / (1 - self::ECC * self::ECC));
        $this->distanceToSun = self::SUNSMAX / $orbitalDistanceFactor;
        $this->sunAngularDiameter = $orbitalDistanceFactor * self::SUNANGSIZ;

        /*
         * Calculation of the Moon's position
         */

        $ml = $this->normalizeAngle(13.1763966 * $dateWithinEpoch + self::MMLONG);    // Moon's mean longitude
        $MM = $this->normalizeAngle($ml - 0.1114041 * $dateWithinEpoch - self::MMLONGP);  // Moon's mean anomaly
        $MN = $this->normalizeAngle(self::MLNODE - 0.0529539 * $dateWithinEpoch);    // Moon's ascending node mean longitude
        $Ev = 1.2739 * sin(deg2rad(2 * ($ml - $eclipticLong) - $MM));  // Evection
        $Ae = 0.1858 * sin(deg2rad($perigeeToEpoch));        // Annual equation
        $A3 = 0.37 * sin(deg2rad($perigeeToEpoch));         // Correction term
        $MmP = $MM + $Ev - $Ae - $A3;         // Corrected anomaly
        $mEc = 6.2886 * sin(deg2rad($MmP));        // Correction for the equation of the centre
        $A4 = 0.214 * sin(deg2rad(2 * $MmP));       // Another correction term
        $lP = $ml + $Ev + $mEc - $Ae + $A4;        // Corrected longitude
        $V = 0.6583 * sin(deg2rad(2 * ($lP - $eclipticLong)));    // Variation
        $lPP = $lP + $V;            // True longitude
        $NP = $MN - 0.16 * sin(deg2rad($perigeeToEpoch));       // Corrected longitude of the node
        $y = sin(deg2rad($lPP - $NP)) * cos(deg2rad(self::MINC));   // Y inclination coordinate
        $x = cos(deg2rad($lPP - $NP));         // X inclination coordinate

        $Lambdamoon = rad2deg(atan2($y, $x)) + $NP;      // Ecliptic longitude
        $BetaM = rad2deg(asin(sin(deg2rad($lPP - $NP)) * sin(deg2rad(self::MINC))));  // Ecliptic latitude

        /* Calculation of the phase of the Moon */
        $MoonAge = $lPP - $eclipticLong;         // Age of the Moon in degrees
        $MoonPhase = (1 - cos(deg2rad($MoonAge))) / 2;     // Phase of the Moon
// Distance of moon from the centre of the Earth
        $MoonDist = (self::MSMAX * (1 - self::MECC * self::MECC)) / (1 + self::MECC * cos(deg2rad($MmP + $mEc)));

        $MoonDFrac = $MoonDist / self::MSMAX;
        $MoonAng = self::MANGSIZ / $MoonDFrac;        // Moon's angular diameter
//$MoonPar = self::MPARALLAX / $MoonDFrac;							// Moon's parallax
// store results
        $this->phase = $this->normalizeAngle($MoonAge) / 360;     // Phase (0 to 1)
        $this->illum = $MoonPhase;          // Illuminated fraction (0 to 1)
        $this->age = self::SYNODIC_MONTH * $this->phase;       // Age of moon (days)
        $this->dist = $MoonDist;          // Distance (kilometres)
        $this->angdia = $MoonAng;          // Angular diameter (degrees)
    }

    /**
     * Normalize angle
     * 
     * @param number $a
     * @return number
     */
    private function normalizeAngle($a)
    {
        return $a - 360 * \floor($a / 360);
    }

    /**
     * Solve Kepler equation
     *  
     * @param number $m
     * @param number $ecc
     * @return number
     */
    private function solveKeplerEquation($m, $ecc)
    {
        $epsilon = \pow(1, -6);

        $e = $m = \deg2rad($m);
        do {
            $delta = $e - $ecc * \sin($e) - $m;
            $e -= $delta / ( 1 - $ecc * \cos($e) );
        } while (\abs($delta) > $epsilon);

        return $e;
    }

    /**
     * Calculate mean new moon
     * 
     * Calculates  time  of  the mean new Moon for a given
     * base date.  This argument K to this function is the
     * precomputed synodic month index, given by:
     *     K = (year - 1900) * 12.3685
     * where year is expressed as a year and fractional year.
     * 
     * @param number $sdate
     * @param number $k
     * @return number
     */
    private function meanphase($sdate, $k)
    {
// Time in Julian centuries from 1900 January 0.5
        $t = ( $sdate - 2415020.0 ) / 36525;
        $t2 = $t * $t;
        $t3 = $t2 * $t;

        $nt1 = 2415020.75933 + self::SYNODIC_MONTH * $k
            + 0.0001178 * $t2
            - 0.000000155 * $t3
            + 0.00033 * sin(deg2rad(166.56 + 132.87 * $t - 0.009173 * $t2));

        return $nt1;
    }

    /**
     * Get true phase time
     * 
     * Given a K value used to determine the mean phase of
     * the new moon, and a phase selector (0.0, 0.25, 0.5,
     * 0.75), obtain the true, corrected phase time.
     * 
     * @param number $k
     * @param number $phase
     * @return boolean|number
     */
    private function truephase($k, $phase)
    {
        $apcor = false;

        $k += $phase;    // Add phase to new moon time
        $t = $k / 1236.85;   // Time in Julian centuries from 1900 January 0.5
        $t2 = $t * $t;    // Square for frequent use
        $t3 = $t2 * $t;    // Cube for frequent use
        $pt = 2415020.75933   // Mean time of phase
            + self::SYNODIC_MONTH * $k
            + 0.0001178 * $t2
            - 0.000000155 * $t3
            + 0.00033 * sin(deg2rad(166.56 + 132.87 * $t - 0.009173 * $t2));

        $m = 359.2242 + 29.10535608 * $k - 0.0000333 * $t2 - 0.00000347 * $t3;   // Sun's mean anomaly
        $mprime = 306.0253 + 385.81691806 * $k + 0.0107306 * $t2 + 0.00001236 * $t3; // Moon's mean anomaly
        $f = 21.2964 + 390.67050646 * $k - 0.0016528 * $t2 - 0.00000239 * $t3;   // Moon's argument of latitude
        if ($phase < 0.01 || abs($phase - 0.5) < 0.01) {
// Corrections for New and Full Moon
            $pt += (0.1734 - 0.000393 * $t) * sin(deg2rad($m))
                + 0.0021 * sin(deg2rad(2 * $m))
                - 0.4068 * sin(deg2rad($mprime))
                + 0.0161 * sin(deg2rad(2 * $mprime))
                - 0.0004 * sin(deg2rad(3 * $mprime))
                + 0.0104 * sin(deg2rad(2 * $f))
                - 0.0051 * sin(deg2rad($m + $mprime))
                - 0.0074 * sin(deg2rad($m - $mprime))
                + 0.0004 * sin(deg2rad(2 * $f + $m))
                - 0.0004 * sin(deg2rad(2 * $f - $m))
                - 0.0006 * sin(deg2rad(2 * $f + $mprime))
                + 0.0010 * sin(deg2rad(2 * $f - $mprime))
                + 0.0005 * sin(deg2rad($m + 2 * $mprime));
            $apcor = true;
        } else if (abs($phase - 0.25) < 0.01 || abs($phase - 0.75) < 0.01) {
            $pt += (0.1721 - 0.0004 * $t) * sin(deg2rad($m))
                + 0.0021 * sin(deg2rad(2 * $m))
                - 0.6280 * sin(deg2rad($mprime))
                + 0.0089 * sin(deg2rad(2 * $mprime))
                - 0.0004 * sin(deg2rad(3 * $mprime))
                + 0.0079 * sin(deg2rad(2 * $f))
                - 0.0119 * sin(deg2rad($m + $mprime))
                - 0.0047 * sin(deg2rad($m - $mprime))
                + 0.0003 * sin(deg2rad(2 * $f + $m))
                - 0.0004 * sin(deg2rad(2 * $f - $m))
                - 0.0006 * sin(deg2rad(2 * $f + $mprime))
                + 0.0021 * sin(deg2rad(2 * $f - $mprime))
                + 0.0003 * sin(deg2rad($m + 2 * $mprime))
                + 0.0004 * sin(deg2rad($m - 2 * $mprime))
                - 0.0003 * sin(deg2rad(2 * $m + $mprime));
            if ($phase < 0.5)  // First quarter correction
                $pt += 0.0028 - 0.0004 * cos(deg2rad($m)) + 0.0003 * cos(deg2rad($mprime));
            else // Last quarter correction
                $pt += -0.0028 + 0.0004 * cos(deg2rad($m)) - 0.0003 * cos(deg2rad($mprime));
            $apcor = true;
        }
        if (!$apcor) // function was called with an invalid phase selector
            return false;

        return $pt;
    }

    /**
     * Find moon phases
     * 
     * Find time of phases of the moon which surround the current date.
     * Five phases are found, starting and
     * ending with the new moons which bound the  current lunation.
     * 
     */
    private function phasehunt()
    {
        $sdate = $this->utcToJulian($this->timestamp);
        $adate = $sdate - 45;
        $ats = $this->timestamp - 86400 * 45;
        $yy = (int) gmdate('Y', $ats);
        $mm = (int) gmdate('n', $ats);

        $k1 = floor(( $yy + ( ( $mm - 1 ) * ( 1 / 12 ) ) - 1900 ) * 12.3685);
        $adate = $nt1 = $this->meanphase($adate, $k1);

        while (true) {
            $adate += self::SYNODIC_MONTH;
            $k2 = $k1 + 1;
            $nt2 = $this->meanphase($adate, $k2);
// if nt2 is close to sdate, then mean phase isn't good enough, we have to be more accurate
            if (abs($nt2 - $sdate) < 0.5)
                $nt2 = $this->truephase($k2, 0.0);
            if ($nt1 <= $sdate && $nt2 > $sdate)
                break;
            $nt1 = $nt2;
            $k1 = $k2;
        }

// results in Julian dates
        $data = array(
            $this->truephase($k1, 0.0),
            $this->truephase($k1, 0.25),
            $this->truephase($k1, 0.5),
            $this->truephase($k1, 0.75),
            $this->truephase($k2, 0.0)
        );

        $this->quarters = array();
        foreach ($data as $v)
            $this->quarters[] = ( $v - 2440587.5 ) * 86400; // convert to UNIX time
    }

    /**
     * Convert timestamp to astronomical Julian time
     * 
     * I.e. Julian date plus day fraction
     * 
     * @param number $ts
     * @return number
     */
    private function utcToJulian($ts)
    {
        return $ts / 86400 + 2440587.5;
    }

    /**
     * Get phase timestamp
     *  
     * @param int $n
     * @return number
     */
    public function getPhaseBeginningDate($n)
    {
        if (\is_null($this->quarters)) {
            $this->phasehunt();
        }

        $dt = \DateTime::createFromFormat('U', (int) $this->quarters[$n]);
        $dt->setTimeZone($this->timezone);

        return $dt;
    }

    /**
     * Get phase [0 to 1]
     * 
     * The terminator phase angle as a fraction of a full circle.
     * Both 0 and 1 correspond to a New Moon, and 0.5 corresponds to a Full Moon.
     * 
     * @return number
     */
    public function phase()
    {
        return $this->phase;
    }

    /**
     * Get illumination fraction [0 to 1]
     * 
     * The illuminated fraction of the Moon (0 = New, 1 = Full).
     * 
     * @return number
     */
    public function illumination()
    {
        return $this->illum;
    }

    /**
     * Get the age of the Moon [d]
     *  
     * @return number
     */
    public function age()
    {
        return $this->age;
    }

    /**
     * Get the distance to the Moon [km]
     * 
     * From the centre of the Earth.
     *  
     * @return number
     */
    public function distance()
    {
        return $this->dist;
    }

    /**
     * Get the Moon angular diameter [rad]
     * 
     * Subtended by the Moon as seen by an observer at the centre of the Earth.
     * 
     * @return number
     */
    public function diameter()
    {
        return $this->angdia;
    }

    /**
     * Get the distance of the Moon to the Sun [km]
     *  
     * @return number
     */
    public function sundistance()
    {
        return $this->distanceToSun;
    }

    /**
     * Get the Sun angular diameter [rad]
     * 
     * Subtended by the Sun as seen by an observer at the centre of the Earth.
     * 
     * @return number
     */
    public function sundiameter()
    {
        return $this->sunAngularDiameter;
    }
    /* Functions for calculating moon rise and moon set */

    /**
     * Calculates moon rise and moon set
     * 
     * Adaptation by the PHP conversion done by Matt Hackmann of an unkown javascript library.
     * 
     * @link http://dxprog.com/entry/calculate-moon-rise-and-set-in-php/
     * @param number $lat
     * @param number $lon
     * @return stdClass->moonrise|moonset
     */
    public function moonrise_and_set($lat, $lon)
    {
        $result = new stdClass;
        $result->moonrise = NULL;
        $result->moonset = NULL;

// Local variables
        $timezone = floor($lon / 15);
        $date = $this->modifiedJulianDate($this->timestamp);
        $date -= ($timezone / 24);
        $latRad = deg2rad($lat);
        $sinho = 0.0023271056;
        $cglat = cos($latRad);
        $sglat = sin($latRad);
        $rise = FALSE;
        $set = FALSE;
        $utrise = 0;
        $utset = 0;

// Loop over hours until the moon rise and set times are found 
        $hour = 1;
        $ym = $this->sinAlt($date, $hour - 1, $lon, $cglat, $sglat) - $sinho;

        while ($hour < 25 and ($rise === FALSE or $set === FALSE)) {
            $yz = $this->sinAlt($date, $hour, $lon, $cglat, $sglat) - $sinho;
            $yp = $this->sinAlt($date, $hour + 1, $lon, $cglat, $sglat) - $sinho;

            list($nz, $z1, $z2, $xe, $ye) = $this->quad($ym, $yz, $yp);

            if ($nz === 1) {
                if ($ym < 0) {
                    $utrise = $hour + $z1;
                    $rise = TRUE;
                } else {
                    $utset = $hour + $z1;
                    $set = TRUE;
                }
            } else if ($nz == 2) {
                if ($ye < 0) {
                    $utrise = $hour + $z2;
                    $utset = $hour + $z1;
                } else {
                    $utrise = $hour + $z1;
                    $utset = $hour + $z2;
                }
            }

            $ym = $yp;
            $hour += 2;
        }

        $utrise = $this->convertDecimalHours($utrise);
        $utset = $this->convertDecimalHours($utset);
        $result->moonrise = $rise ? strtotime("Midnight +{$utrise['hrs']} hours +{$utrise['min']} minutes",
                $this->timestamp) : strtotime('Midnight +1 day', $this->timestamp);
        $result->moonset = $set ? strtotime("Midnight +{$utset['hrs']} hours +{$utset['min']} minutes",
                $this->timestamp) : strtotime('Midnight +1 day', $this->timestamp);

        return $result;
    }

    /**
     * Caculates the sine of the moon altitude
     * 
     * @param number $mjd (Julian date) [s]
     * @param number $hour (Hour) [h]
     * @param number $glon (Longitude) [deg]
     * @param number $cglat (Cosine of the latitude)
     * @param number $sglat (Sine of the latitude)
     * @return number
     */
    private function sinAlt($mjd, $hour, $glon, $cglat, $sglat)
    {
        $mjd += $hour / 24;
        $t = ($mjd - 51544.5) / 36525;
        $objpos = $this->minimoon($t);

        $ra = $objpos[1];
        $dec = $objpos[0];
        $decRad = deg2rad($dec);
        $tau = 15 * ($this->lmst($mjd, $glon) - $ra);

        return $sglat * sin($decRad) + $cglat * cos($decRad) * cos(deg2rad($tau));
    }

    /**
     * Takes t and returns the geocentric ra and dec in an array
     * 
     * Claimed good to 5' (angle) in ra and 1' in dec
     * Tallies with another approximate method and with ICE for a couple of dates
     * 
     * @param number $t
     * @return array($dec, $ra)
     */
    private function minimoon($t)
    {
        $p2 = 6.283185307;
        $arc = 206264.8062;
        $coseps = 0.91748;
        $sineps = 0.39778;

        $lo = $this->frac(0.606433 + 1336.855225 * $t);
        $l = $p2 * $this->frac(0.374897 + 1325.552410 * $t);
        $l2 = $l * 2;
        $ls = $p2 * $this->frac(0.993133 + 99.997361 * $t);
        $d = $p2 * $this->frac(0.827361 + 1236.853086 * $t);
        $d2 = $d * 2;
        $f = $p2 * $this->frac(0.259086 + 1342.227825 * $t);
        $f2 = $f * 2;

        $sinls = sin($ls);
        $sinf2 = sin($f2);

        $dl = 22640 * sin($l);
        $dl += -4586 * sin($l - $d2);
        $dl += 2370 * sin($d2);
        $dl += 769 * sin($l2);
        $dl += -668 * $sinls;
        $dl += -412 * $sinf2;
        $dl += -212 * sin($l2 - $d2);
        $dl += -206 * sin($l + $ls - $d2);
        $dl += 192 * sin($l + $d2);
        $dl += -165 * sin($ls - $d2);
        $dl += -125 * sin($d);
        $dl += -110 * sin($l + $ls);
        $dl += 148 * sin($l - $ls);
        $dl += -55 * sin($f2 - $d2);

        $s = $f + ($dl + 412 * $sinf2 + 541 * $sinls) / $arc;
        $h = $f - $d2;
        $n = -526 * sin($h);
        $n += 44 * sin($l + $h);
        $n += -31 * sin(-$l + $h);
        $n += -23 * sin($ls + $h);
        $n += 11 * sin(-$ls + $h);
        $n += -25 * sin(-$l2 + $f);
        $n += 21 * sin(-$l + $f);

        $L_moon = $p2 * $this->frac($lo + $dl / 1296000);
        $B_moon = (18520.0 * sin($s) + $n) / $arc;

        $cb = cos($B_moon);
        $x = $cb * cos($L_moon);
        $v = $cb * sin($L_moon);
        $w = sin($B_moon);
        $y = $coseps * $v - $sineps * $w;
        $z = $sineps * $v + $coseps * $w;
        $rho = sqrt(1 - $z * $z);
        $dec = (360 / $p2) * atan($z / $rho);
        $ra = (48 / $p2) * atan($y / ($x + $rho));
        $ra = $ra < 0 ? $ra + 24 : $ra;

        return array($dec, $ra);
    }

    /**
     * Returns the fractional part of x
     *
     * @param number $x
     * @return number
     */
    private function frac($x)
    {
        $x -= (int) $x;
        return $x < 0 ? $x + 1 : $x;
    }

    /**
     * Calculates LMST?
     * 
     * @param number $mjd (Julian time) [s]
     * @param number $glon (Longitude) [deg]
     * @return number
     */
    private function lmst($mjd, $glon)
    {
        $d = $mjd - 51544.5;
        $t = $d / 36525;
        $lst = $this->degRange(280.46061839 + 360.98564736629 * $d + 0.000387933 * $t * $t - $t * $t * $t / 38710000);
        return $lst / 15 + $glon / 15;
    }

    /**
     * Returns an angle in degrees between 0 and 360
     * 
     * @param number $x
     * @return number
     */
    private function degRange($x)
    {
        $b = $x / 360;
        $a = 360 * ($b - (int) $b);
        $retVal = $a < 0 ? $a + 360 : $a;
        return $retVal;
    }

    /**
     * Solve quadratic equation
     * 
     * Finds the parabola throuh the three points (-1,ym), (0,yz), (1, yp)
     * and returns the coordinates of the max/min (if any) xe, ye
     * the values of x where the parabola crosses zero (roots of the self::quadratic)
     * and the number of roots (0, 1 or 2) within the interval [-1, 1]
     *
     * Well, this routine is producing sensible answers
     *
     * Results passed as array [nz, z1, z2, xe, ye]
     * 
     * @param number $ym
     * @param number $yz
     * @param number $yp
     * @return array
     */
    private function quad($ym, $yz, $yp)
    {

        $nz = $z1 = $z2 = 0;
        $a = 0.5 * ($ym + $yp) - $yz;
        $b = 0.5 * ($yp - $ym);
        $c = $yz;
        $xe = -$b / (2 * $a);
        $ye = ($a * $xe + $b) * $xe + $c;
        $dis = $b * $b - 4 * $a * $c;
        if ($dis > 0) {
            $dx = 0.5 * sqrt($dis) / abs($a);
            $z1 = $xe - $dx;
            $z2 = $xe + $dx;
            $nz = abs($z1) < 1 ? $nz + 1 : $nz;
            $nz = abs($z2) < 1 ? $nz + 1 : $nz;
            $z1 = $z1 < -1 ? $z2 : $z1;
        }

        return array($nz, $z1, $z2, $xe, $ye);
    }

    /**
     * Calculate modified julian date
     * 
     * Takes the day, month, year and hours in the day and returns the
     * modified julian day number defined as mjd = jd - 2400000.5
     * checked OK for Greg era dates - 26th Dec 02
     * 
     * @param $timestamp
     * @return $juliandate
     */
    private function modifiedJulianDate($timestamp)
    {
        $day = gmdate('d', $timestamp);
        $month = gmdate('m', $timestamp);
        $year = gmdate('Y', $timestamp);

        if ($month <= 2) {
            $month += 12;
            $year--;
        }

        $a = 10000 * $year + 100 * $month + $day;
        $b = 0;
        if ($a <= 15821004.1) {
            $b = -2 * (int) (($year + 4716) / 4) - 1179;
        } else {
            $b = (int) ($year / 400) - (int) ($year / 100) + (int) ($year / 4);
        }

        $a = 365 * $year - 679004;
        return $a + $b + (int) (30.6001 * ($month + 1)) + $day;
    }

    /**
     * Convert decimal hours to hours and minutes
     * 
     * @param float $hours
     * @return array('hrs'|'min')
     */
    private function convertDecimalHours($hours)
    {
        $hrs = (int) ($hours * 60 + 0.5) / 60.0;
        $h = (int) ($hrs);
        $m = (int) (60 * ($hrs - $h) + 0.5);
        return array('hrs' => $h, 'min' => $m);
    }
}
