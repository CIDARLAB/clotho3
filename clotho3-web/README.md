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

- Fixes
	- searchbar posts all say()s

- Editor
	- editor parent, with children for specifics (function, schema, trail, etc.)
	- make trail editor template
	- better dynamic editor    
	- function editor tests reset
	- gradeQuiz
    		- JSON stringify results
    		- use for function tests

- Trail Content
	- easier to edit, organize

- App organization
    - get rid of Application prefixing of angular modules
    - break up ui components
    - add grunt, bower, yeoman support




- searchbar updates
    - call say() for run(), other commands
        - check isString in say and wrap

- clotho.run directive to handle module namespacing, or use submit

- DIY construction file


DNA functions
    - ligation alignment, HTML
    - handle overhangs in anneal


- trail: unread check for each Page, change style, affects dialogs

- $focus.typeOutSearch to accept array, use in Ligation

- construction file -> protocol




-------------------------


- $focus
    - $dialog with custom top spacing
    - pulse submit button
    - move search typeout to search service?
    - fn: show activity log (use searchbar service)
    - type series of commands with break


- custom person object, add anderson lab people
    - start testing andersonLab site with that (ninaEmami ? )


- new schemas
    - Usage
        - e.g. for function, construction file, etc.
    - Construction File
        - Schemas for steps (digest, pcr, ligation, and all subclasses)


- UI services
    - $focus
        - get popover to display by default
    - help tips service (as popoverService)
        - given object of tips and selectors, show appropriately
        - be able to step through
    - nice dropdown for filtering
          http://jsfiddle.net/TahmidTanzim/N9Vqk/



- $timeout API requests e.g. get, query, after interval


- upgrade to $modal from angular ui


- finish schema editor -- verify with stephanie
    - handle parameterized fields (e.g. array | module)
    - methods section
        - custom template to show id?
        - way of adding tags (http://jsfiddle.net/joshdmiller/hAz5A/)
    - selecting parent should show parent fields (small and disabled)
    - collect after issue query, save, update editor


- DNA Functions
    - decide data structure
    - DNA
        - required upstream
            - species specific
                - codon optimize
                - silent sites
        - parse: gb, fasta, embl
        - BLAST breakout
    - Digest
        - import NEB enzymes, use same format
            - broaden methylation categories: dam, dcm, cpg
            - import
                - formats: ftp://ftp.neb.com/pub/rebase/REBASE.DOC
                - GCG format: ftp://ftp.neb.com/pub/rebase/gcg.txt
                - also: star, homing, etc. etc. etc. -- write parsers
        - dep on DNA service new functions
            - swapout Digest site (need to know ORF) -- not demo
        - import strider format (REBASE) / NEB format
    - PCR
        - other pcr types


- anderson lab site
    - slideshow directive
    - people page with ng-repeat



- Command Bar rewrite
    - note -- inherit some UI stuff from menuController so not in searchbar
    - use typeahead with custom template??
        - want to reduce bootstrap reliance, but well-written. maybe move to core
    - add login functionality (not showing for demo)
    - pass in options
    - rename CSS classes
    - avoid ui-directives.js etc.
    - reduce template use


-typeahead
    - Clotho.query() promise as async match -- get working when not in controller
    - custom template: http://stackoverflow.com/questions/18245834/bootstrap-ui-typeahead-display-more-than-one-property-in-results-list/18251561#18251561



- API
    - collect queries
    - verify collector updates on set() and destroy() and get()
    - Socket events in own layer (outside PubSub)
    - add-in $dialog for alert() and edit() but don't make dependency
        - remove from anderson_lab dependencies


- Trails
    - youtube directive
        - extract param for autoplay in attrs
        - slideout attribute (show only when request, don't autoplay)
        - listen for changes to id, recompile
    - tool -- how to load?
        - clotho-run is easy...
        - widgets?
    - demo event at time in video (e.g. type something via jQuery)
        -- also demonstrate can use angular $scope events (like next() and will compile)
    - Trail schema writeup / Trail
    - Trail editor -- req. schema finalized
    - prettyprint directive
    - Quizzes
        - answer -> color side nav
        - templates
            - write multi-tf (need directive?)
            - multipic : active class



- Mac Housekeeping
    - folder actions to push changes to models folder to testData folder
    - tie github pulls to original repositories and keep in sync
    - package clotho api stuff (foundation folder and application.js)
        - grunt / yeoman --- javascript project builder


-


- editor
    - break up controllers in generic-editor-directive
        - reset tests on changing id, editMode, etc. -- do once separate controllers
    - start edit mode on createSchema
    - replace internal HTML, not whole thing (closeable directive)
    - use "ID" not "UUID"
    - use accordion for tests div
    - map BSON -> HTML5 data types
    - query for test args
        - show name, not uuid
        - Popover with whole JSON shown
    - use form controller for validation
        - currently requires nested forms: https://github.com/angular/angular.js/issues/1404
    - new fields
        - multiple select
        - select
        - radios
        - file drops?
    - add novalidate via directive
        - need to test novalidate="true"
        - custom validation using ngForm or ngModel.$setValidity
    - function editor
        - integrate codemirror, tie to language selection
    - deep linking using search and $routeProvider.when( ... reloadOnSearch = false)


- UI Directives
    - popover with HTML (see angularStrap, or wait for bootstrap UI)
        - https://github.com/angular-ui/bootstrap/issues/220
        - hide on mouseleave when mouseenter trigger
        - angular Strap: https://raw.github.com/mgcrea/angular-strap/v0.7.5/dist/angular-strap.js
    - rewrite contenteditable
        - make play nicely (e.g. with digest-highlight)
    - DRAGGABLE
        - http://www.bennadel.com/blog/2446-Using-Controllers-In-Directives-In-AngularJS.htm
        - https://github.com/fisshy/Angular-drag-drop
        - https://github.com/fisshy/Angular-drag-drop/blob/master/javascript/dragandrop.js
        - http://jsfiddle.net/ADukg/2516/
        - need to handle drop + other events


- plasmid editors to check out
    - benchling
    - sparkDNA



- caret movement directives : https://github.com/DrPheltRight/jquery-caret
    - see also contenteditable directive - have one working
    - get current position and save as option



- bootstrapping apps without ng-view (just use ng-include or templates)
    - use ng 1.2.0 - routes not included by default
    - register app id with module


- break up UI directives, filters, etc. into modules



- (cool idea // advantage?) plasmid editor -- add locations parser to ngModel pipeline, store locations inside of ngModel (on the object), push them back in using a formatter



- $keypress.suspend() e.g. in trail mode when using Command Bar -- need to be namespaced


- favorite() API method - and unfavorite() and dislike()



UI STUFF
- services to write / borrow
    - progressBar (bootstrap)
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




- Refine PubSub / communication
    - REMOVE SOCKET LISTENERS (only listen to PubSub)
    - refine pubsub listeners object
    - PubSub clear needs to clean 'references' object
    - references
        - http://closure-library.googlecode.com/svn/docs/closure_goog_pubsub_pubsub.js.source.html
        - http://www.gridlinked.info/angularJS/modules/MessagingServices.js (requires jQuery)



- Routing
    - test asynchronously adding routes via controller
        - reference
           - http://stackoverflow.com/questions/13153121/how-to-defer-routes-definition-in-angular-js
           - https://groups.google.com/forum/#!msg/angular/mrcy_2BZavQ/Mqte8AvEh0QJ
           - http://stackoverflow.com/questions/13681116/angularjs-dynamic-routing





---------- maybes --------------


- (?) write ClothoProvider - include in module config
    - look for other advantages / utility



- decorators? i.e. for $http


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