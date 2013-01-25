<?php

namespace Swis\Bundle\MoonPhaseBundle\Service;

use Swis\Bundle\MoonPhaseBundle\Model\MoonPhase;

/**
 * Moon phase calculation class.
 * 
 * Originally adapted for PHP from Moontool for Windows (http://www.fourmilab.ch/moontoolw/) by Samir Shah.
 * Transformed and bundled for Symfony by Sascha Weyers.
 *
 * @link https://github.com/swis/MoonPhaseBundle
 * @link https://github.com/solarissmoke/php-moon-phase
 * @author Samir Shah (http://rayofsolaris.net)
 * @author Sascha Weyers
 * @license MIT
 */
class MoonPhaseCalculator
{

    private $synodicMonth = null;

    public function __construct()
    {
        $this->synodicMonth = new \DateInterval('P29DT12H44M3S');
    }

    public function getLastNewMoon(\DateTime $dt = null)
    {
        $moonPhase = new MoonPhase($dt);
        return $moonPhase->getPhaseBeginningDate(0);
    }

    public function getNextNewMoon(\DateTime $dt = null)
    {
        $moonPhase = new MoonPhase($dt);
        return $moonPhase->getPhaseBeginningDate(4);
    }

    public function getLastFullMoon(\DateTime $dt = null)
    {
        $moonPhase = new MoonPhase($dt);
        $fullMoon = $moonPhase->getPhaseBeginningDate(2);
        
        if ($fullMoon > $dt) {
            $fullMoon->sub($this->synodicMonth);
        }
        
        return $fullMoon;
    }
    
    public function getNextFullMoon(\DateTime $dt = null)
    {
        $moonPhase = new MoonPhase($dt);
        $fullMoon = $moonPhase->getPhaseBeginningDate(2);
        
        if ($fullMoon < $dt) {
            $fullMoon->add($this->synodicMonth);
        }
        
        return $fullMoon;
    }
}
