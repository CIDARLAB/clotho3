# Clotho 3.0 Sharable Editor - AngularJS

## Overview

Clotho 3.0 Front-end Playground: Functionalities for:
- Sharable Editor
- Trails Browser
- Chat (socket implementation)

## To Run

- set current directory to this one, then run `node scripts/web-server-socket.js` and go to localhost:8080/app/index.html
- note: socket.io must be installed locally. See http://socket.io/#how-to-use

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

#### Notes
- https://github.com/angular/angular.js/wiki/The-Nuances-of-Scope-Prototypal-Inheritance

### Personal notes

#### Compatibility
- avoid directive element names (type E) for IE compatibility, or use document.createElement (i.e. for non-native elements like <widget></widget>)
    - see AngularUI ie-shiv


## TODOs, Future references, etc.

### Important Fixes + Changes



### Current Tasks

- grunt / yeoman --- javascript project builder
    - package clotho api stuff (foundation folder and application.js)`


- Server communication
    - create custom layer for Socket Events (using PubSub alternative implementation)
    - basic node versions of:
        - query
        - run
        - create
        - revert


- DNA Functions
    - decide data structure
    - DNA
        - determine ORF
        - parse: gb, fasta, embl
        - silent sites
        - BLAST
    - Digest
        - swapout Digest site
        - import strider format (REBASE) / NEB format
        - Dam/Dcm methylation
    - PCR
        - other pcr types


- create help tips - use popovers? - structure JSON


- $clotho.api
    - ensure only instantiated once


- editor
    - use form controller for validation
        - currently requires nested forms: https://github.com/angular/angular.js/issues/1404
    - add novalidate via directive
    - new fields
        - select
        - radios
        - file drops?

    - function editor
        - integrate codemirror, tie to language selection

    - schema editor -- wait on stephanie



- UI Directives
    - typeahead with async content
        - http://stackoverflow.com/questions/15930339/how-to-tie-angular-uis-typeahead-with-a-server-via-http-for-server-side-optimi
        - http://plnkr.co/edit/eGG9Kj?p=preview
    - rewrite contenteditable
        - make play nicely (e.g. with digest-highlight)
    - dropdown in Angular --  see angular-strap
    - popover with HTML -- test -- see angular-strap
        - hide on mouseleave when mouseenter trigger
    - angular Strap: https://raw.github.com/mgcrea/angular-strap/v0.7.5/dist/angular-strap.js



- plasmid editors to check out
    - benchling
    - sparkDNA



- Application Mixin promises
    - terminal controller -- required as dependency
    - youtube API @ Trails -- not required immediately


- caret movement directives : https://github.com/DrPheltRight/jquery-caret
    - see also contenteditable directive - have one working
    - get current position and save as option



- bootstrapping apps without ng-view (just use ng-include)


- create wrappers for directives, filters -> expose DNA functions accordingly



- break up UI directives, filters, etc. into modules




- streamline display and show() + adding components
    - get bootstrapping to work



- migrate to $log from console.log



- draggable elements (within angular -- also sortable)
    - https://github.com/fisshy/Angular-drag-drop/blob/master/javascript/dragandrop.js
    - http://jsfiddle.net/ADukg/2516/



- (advantage?) plasmid editor -- add locations parser to ngModel pipeline, store locations outside of ngModel, push them back in using a formatter



- $keypress.suspend() e.g. in trail mode when using search bar -- need to be namespaced


- move searchbar to typeahead


- trails functionality
    - youtube events - pass as param to function
    - keybindings for advancing etc.
    - templates
        - types
            - sequence
            - assertion / reason
            - hot spot (choose the wrong thing) /// paragraph with multi TFs
            - multipic - add active class
            - drag-drop - make draggable
            - ranking - tie to model
    - UI directives
        - draggable
            - see UI-sortable
            - drag-to for drag-drop
    - answering quiz colors + disables side-nav


- Search Bar
    - reset autocomplete detail when changing
    - animate showing + hiding -- use $animate
    - autocomplete
        - use cache for results?
        - $timeout sending to 50ms?
    - mouseenter and leave events -- maybe size list items differently
        - undetail propagate to children? why?
        - see: http://plnkr.co/edit/QSgRu4Sr8TSqR86eAfvq?p=preview
    - future: rewrite
        - easy to transclude, with customizable / default options



- favorite() API method - and unfavorite() and dislike()


- clickOutside handlers - better logic esp. for ID


- Dynamic Adds
    - Reference
        - watch : https://github.com/mhevery/angular.js/blob/dte/src/auto/injector.js
        - controllers : https://github.com/matys84pl/angularjs-requirejs-lazy-controllers
    - load and apply css (via jquery? - whatever use for getting scripts)
        - test downloading css and see if applies
    - simple Clotho.show()
        - get one master object for controller
        - template tied to that controller




UI STUFF
- sharable modal --- make it actually share
- services to write / borrow
    - dropdown (bootstrap)
        - use for searchbar help
    - progressBar (bootstrap)
    - scrollTo (probably easy to write) / scrollfix (see angular ui utils)
- ANIMATION
    - see CSS framework https://github.com/daneden/animate.css


- Clotho.show()
    - store model and template name somewhere, access in dynamic_template for show()
        - after download JS: callback or something after controller registered?
        - Register as controller object etc. rather than global function
            - EXAMPLE SEE http://stackoverflow.com/questions/15250644/angularjs-loading-a-controller-dynamically
        - Caching templates
            - remove and reset if exists
                - http://deansofer.com/posts/view/14/AngularJs-Tips-and-Tricks-UPDATED#refresh-partials


    - show() in modal box as option --- make this alert?
        - simple example: http://jsfiddle.net/jjzabkar/MWAek/1/
        - see angular UI / UI-bootstrap modal directive
            - http://angular-ui.github.io/bootstrap/




- TESTING
    - write unit tests, etc.




- Clotho API
    - once REST works, write synchronous Clotho.get()
    - combine watch and watch2 - check for function vs. obj and field
        - try using angular.extend()



- remove listeners
    - see $scope.$watch:
        - deregistration returned as function to invoke later
        - http://www.bennadel.com/blog/2480-Unbinding-watch-Listeners-In-AngularJS.htm
    - reg app ID with each app
    - make universal --> $scope.$on('$destroy', Clotho.silence($scope.$id))
            -i.e. (?) extend $rootScope.$destroy() if possible



- Refine PubSub / communication
    - REMOVE SOCKET LISTENERS (only listen to PubSub)
    - refine pubsub listeners object
    - PubSub clear needs to clean 'references' object
    - references
        - http://closure-library.googlecode.com/svn/docs/closure_goog_pubsub_pubsub.js.source.html
        - http://www.gridlinked.info/angularJS/modules/MessagingServices.js (requires jQuery)



- Routing
    - app with no routes (i.e. for widgets) ???
    - test asynchronously adding routes via controller
        - reference
           - http://stackoverflow.com/questions/13153121/how-to-defer-routes-definition-in-angular-js
           - https://groups.google.com/forum/#!msg/angular/mrcy_2BZavQ/Mqte8AvEh0QJ
           - http://stackoverflow.com/questions/13681116/angularjs-dynamic-routing



- modal & $dialog
    - create mode for Clotho API
        - assign uuid to dialog (should come from server?)
    - fix $animation (on way out) (or wait for project to do it)




- deep-linking for search bar and editor
    - http://stackoverflow.com/questions/14974271/can-you-change-a-path-without-reloading-the-controller-in-angularjs





-Editor
    - tied to ngForm
    - form coloring - check for both: required (red), and validation (yellow)
        - class assignment not working as expected right now... waiting for next angular steady version
        - http://stackoverflow.com/questions/7792652/what-is-the-best-way-to-conditionally-apply-a-class-with-angularjs/8309832#8309832

    - file drop as input type

    - SCHEMA
        - add hierarchies in schema?
        - permissions - pass in sharable data, not schema

    - REGEXP FOR VALIDATION e.g. state
        - use ng-pattern: use regexp e.g. $scope.state =  /^\w\w$/





---------- maybes --------------


- (?) write ClothoProvider - include in module config
    - look for other advantages / utility



- decorators? i.e. for $http



- schema builder


------ experiments ------

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

#### ActiveMQ & node.js

- http://stackoverflow.com/questions/4700935/what-are-good-message-queue-options-for-nodejs

#### PubNub

- http://cdn.pubnub.com/pubnub-3.1.js
- http://jsfiddle.net/bv5Kq/12/

#### Authentication

- http://www.espeo.pl/2012/02/26/authentication-in-angularjs-application

#### Widgets

- http://jsfiddle.net/simpulton/VJ94U/
- http://onehungrymind.com/angularjs-sticky-notes-pt-1-architecture/
- Simple Draggable Directive: http://docs.angularjs.org/guide/compiler
- Draggable things in Angular (with shared controller w/o require): http://www.bennadel.com/blog/2446-Using-Controllers-In-Directives-In-AngularJS.htm

#### Manual Bootstrapping / Multiple Apps per page

- http://docs.angularjs.org/api/angular.bootstrap
- https://groups.google.com/forum/#!msg/angular/lhbrIG5aBX4/4hYnzq2eGZwJ