<?php

namespace Swis\Bundle\MoonPhaseBundle\Tests\Service;

use Swis\Bundle\MoonPhaseBundle\Service\MoonPhaseCalculator;
use Symfony\Bundle\FrameworkBundle\Test\WebTestCase;

class MoonPhaseCalculatorTest extends WebTestCase
{

    public function testMoonPhaseCalculation()
    {
        $calc = new MoonPhaseCalculator();

        $dt = \DateTime::createFromFormat('YmdHis', '20130126000000');
        $dt->setTimeZone(new \DateTimeZone('UTC'));

        $this->assertEquals(201301111945, $calc->getLastNewMoon($dt)->format('YmdHi'));
        $this->assertEquals(201302100722, $calc->getNextNewMoon($dt)->format('YmdHi'));

        $this->assertEquals(201212281556, $calc->getLastFullMoon($dt)->format('YmdHi'));
        $this->assertEquals(201301270440, $calc->getNextFullMoon($dt)->format('YmdHi'));

        $dt = \DateTime::createFromFormat('YmdHis', '20120802000000');
        $dt->setTimeZone(new \DateTimeZone('UTC'));

        $this->assertEquals(201207190424, $calc->getLastNewMoon($dt)->format('YmdHi'));
        $this->assertEquals(201208171554, $calc->getNextNewMoon($dt)->format('YmdHi'));

        $this->assertEquals(201207031443, $calc->getLastFullMoon($dt)->format('YmdHi'));
        $this->assertEquals(201208020327, $calc->getNextFullMoon($dt)->format('YmdHi'));
    }
    
    public function test()
    {
        $calc = new MoonPhaseCalculator();
        
        $dt = \DateTime::createFromFormat('YmdHis', '20120802000000');
        $dt->setTimeZone(new \DateTimeZone('UTC'));

        $fúllMoonOne = $calc->getLastNewMoon($dt);
        $fullMoonTwo = $calc->getNextNewMoon($fúllMoonOne);
        
        $this->assertEquals($fúllMoonOne->format('YmdHi'), $fullMoonTwo->format('YmdHi'));
    }
}
