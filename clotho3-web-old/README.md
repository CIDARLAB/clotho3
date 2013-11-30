# Clotho 3.0 Sharable Editor - AngularJS

## Overview

Clotho 3.0 Front-end Playground: Functionalities for:
- Sharable Editor
- Trails Browser
- Chat (socket implementation) 

## Prerequisites

### Git

### Node.js
- web-server currently runs on node.

### socket.io
- Must be installed for Socket communication. See: http://socket.io/#how-to-use

### Java
- http://www.java.com




## References

### From Angular

- DOCS
- API

#### Angular Notes
- https://github.com/angular/angular.js/wiki/The-Nuances-of-Scope-Prototypal-Inheritance

### Personal notes


## TODOs, Future references, etc.

### Important Fixes + Changes



### Current Tasks

 --- SEE GITHUB ISSUE TRACKING ----



LOW PRIORITY

- UI services
    - nice dropdown for filtering
          http://jsfiddle.net/TahmidTanzim/N9Vqk/


- upgrade to $modal from angular ui


- UI Directives
    - popover with HTML (see angularStrap, or wait for bootstrap UI)
        - https://github.com/angular-ui/bootstrap/issues/220
        - hide on mouseleave when mouseenter trigger
        - angular Strap: https://raw.github.com/mgcrea/angular-strap/v0.7.5/dist/angular-strap.js
    - rewrite contenteditable
        - make play nicely (e.g. with digest-highlight)


- Clotho API
    - once REST works, write synchronous Clotho.get()

------ experiments ------

- (cool idea // advantage?) plasmid editor -- add locations parser to ngModel pipeline, store locations inside of ngModel (on the object), push them back in using a formatter

- filters that perform an action / broadcast events / etc.



#### Research

- re-read: http://deansofer.com/posts/view/14/AngularJs-Tips-and-Tricks-UPDATED#refresh-partials
- currying? angular.bind? use cases?

### Future tasks

- D3.js to make simple plasmids

### TODOs by category


#### General UI

- alerts via Bootstrap (and modals?)

#### General Backend


#### Forward looking


### At Deployment

#### General

- minify and concatenate JS
- non-cacheable template
- cache with versions: js, css, views, code
- sanitize (escape) HTML
- IE <7 compatibility
    - document.createElement(<custom tag names>)
    - JSON polyfill

#### AJAX Crawlable

- https://developers.google.com/webmasters/ajax-crawling/docs/specification

## Things / Apps to check out

#### General Examples

- https://github.com/angular/angular.js/wiki/JsFiddle-Examples
- http://builtwith.angularjs.org/

#### Authentication

- http://www.espeo.pl/2012/02/26/authentication-in-angularjs-application

