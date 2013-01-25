MoonPhaseBundle
===============
[![Build Status](https://travis-ci.org/swis/MoonPhaseBundle.png)](https://travis-ci.org/swis/MoonPhaseBundle)

This Bundle is designed to add some moon phase calculation stuff to your Symfony2 project.

## Installation

### Installation by Composer

If you use composer, add the bundle as dependency to the composer.json of your application

#### for Symfony 2.x

    "require": {
        ...
        "swis/moonphase-bundle": "dev-master"
        ...
    },

### Installation via git submodule

```bash
    git submodule add git://github.com/swis/MoonPhaseBundle.git vendor/bundles/Swis/Bundle/MoonPhaseBundle
```

## Instantiate Bundle in your kernel init file

```php
// app/AppKernel.php
<?php
    // ...
    public function registerBundles()
    {
        $bundles = array(
            // ...
            new Swis\Bundle\MoonPhaseBundle\SwisMoonPhaseBundleBundle(),
        );
    }
```

## Base configuration

Actually, no configuration needed.

## Usage

In your controller, use

```php
<?php
    $this->get('swis_moon_calculator')->getLastNewMoon($dt);
```
