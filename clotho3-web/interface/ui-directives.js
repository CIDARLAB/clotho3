'use strict';

Application.Interface.directive.code = function() {
    return { restrict:'E', terminal: true };
};

/***********************
 CLICKS / COMMON BROWSER EVENTS
 ***********************/

/**
 * Taken from Angular UI
 *
 * General-purpose Event binding. Bind any event not natively supported by Angular
 * Pass an object with keynames for events to ui-event
 * Allows $event object and $params object to be passed
 *
 * @example <input ui-event="{ focus : 'counter++', blur : 'someCallback()' }">
 * @example <input ui-event="{ myCustomEvent : 'myEventHandler($event, $params)'}">
 *
 * @param ui-event {string|object literal} The event to bind to as a string or a hash of events with their callbacks
 */
Application.Interface.directive('uiEvent', ['$parse',
    function ($parse) {
        return function ($scope, elm, attrs) {
            var events = $scope.$eval(attrs.uiEvent);
            angular.forEach(events, function (uiEvent, eventName) {
                var fn = $parse(uiEvent);
                elm.bind(eventName, function (evt) {
                    var params = Array.prototype.slice.call(arguments);
                    //Take out first paramater (event object);
                    params = params.splice(1);
                    fn($scope, {$event: evt, $params: params});
                    if (!$scope.$$phase) {
                        $scope.$apply();
                    }
                });
            });
        };
    }]);

// note - requires jQuery
// future - to make more universal, see:
// https://raw.github.com/cowboy/jquery-outside-events/v1.1/jquery.ba-outside-events.js
// angular way of creating ng-directives
// note - Requires manual activiation to avoid too many $digests()
//      - activate with $scope.$broadcast('$event:$active')
Application.Interface.directive('clickOutside', ['$document', '$parse', function($document, $parse) {
    return function(scope, element, attr, ctrl) {
        var handler = function(event) {
            //todo - rewrite without jQuery.has()
            if (element.has(event.target).length == 0)
                scope.$safeApply(
                    $parse(attr['clickOutside'])(scope, {$event: event})
                );
        };

        //kill on scope desctruction
        scope.$on('$destroy', function() {
            $document.unbind('click', handler);
        });

        //custom events
        scope.$on("clickOutside:$active", function(event, id) {
            if (scope.$id == id) {
                event.preventDefault();
                $document.bind('click', handler);
            }
        });

        scope.$on("clickOutside:$inactive", function(event, id) {
            if (scope.$id == id) {
                event.preventDefault();
                $document.unbind('click', handler);
            }
        });
    }
}]);

//todo - test -- prevent on mouseleave child element -- prevent propagation
Application.Interface.directive('myMouseleave', ['$parse', function($parse) {
    return function(scope, element, attr) {
        element.bind('mouseleave', function(event) {
            event.stopPropagation();
            $parse(attr['myMouseleave'])(scope, {$event: event});
        });
    };
}]);

//future - remove polyfill once available-- currently not present in angular (1.1.5)
Application.Interface.directive('ngFocus', ['$parse', function($parse) {
    return function(scope, element, attr) {
        var fn = $parse(attr['ngFocus']);
        element.bind('focus', function(event) {
            scope.$apply(function() {
                fn(scope, {$event:event});
            });
        });
    }
}]);

//future - remove polyfill once available-- currently not present in angular (1.1.5)
Application.Interface.directive('ngBlur', ['$parse', function($parse) {
    return function(scope, element, attr) {
        var fn = $parse(attr['ngBlur']);
        element.bind('blur', function(event) {
            scope.$apply(function() {
                fn(scope, {$event:event});
            });
        });
    }
}]);

/***********************
 * HTML5 EXTENSIONS
***********************/

// todo - handle val() if present, default to text()

Application.Interface.directive('contenteditable', [function() {

    //moves cursor to end of contenteditable element -- not textarea (but those should be automatic)
    //kinda a hack, but gives same behavior as input
    //todo - move to $caret service
    //future - get cursor location, and reset to there rather than end
    //note - expects DOM element, not jQuery
    function setEndOfContenteditable(contentEditableElement)
    {
        var range,selection;
        if(document.createRange)//Firefox, Chrome, Opera, Safari, IE 9+
        {
            range = document.createRange();
            range.selectNodeContents(contentEditableElement);
            range.collapse(false);
            selection = window.getSelection();
            selection.removeAllRanges();
            selection.addRange(range);
        }
        else if(document.selection)//IE 8 and lower
        {
            range = document.body.createTextRange();
            range.moveToElementText(contentEditableElement);
            range.collapse(false);
            range.select();
        }
    }

    return {
        require: '?ngModel',
        link: function(scope, elm, attrs, ctrl) {
            if(!ctrl) return;

            // view -> model
            elm.bind('keyup', function(event) {
                var key = event.keyCode;
                if (key === 91 || (15 < key && key < 19) || (37 <= key && key <= 40)) return;

                var value = elm.text();
                console.log(value);

                //testing
                //console.log("EDITABLE", ctrl.$viewValue, ctrl.$modelValue, value, ctrl.$parsers, ctrl);

                //todo - use parser here, convert &nbsp; to ' '
                //trim
                value = (angular.isString(value) ? value.replace(/^\s*/, '').replace(/\s*$/, '') : value);

                if (ctrl.$viewValue != value) {
                    scope.$apply(function() {
                        ctrl.$setViewValue(value);
                        //testing
                        //console.log("EDITABLE 2", ctrl.$viewValue, ctrl.$modelValue, value, ctrl);
                        setEndOfContenteditable(elm[0]);
                    });
                }
            });

            // model -> view
            ctrl.$render = function() {
                //if (ctrl.$viewValue != elm.text())
                    //elm.text(isEmpty(ctrl.$viewValue) ? '' : ctrl.$viewValue);
                    elm.text(ctrl.$viewValue);
            };

            // load init value from DOM
            //elm.text(ctrl.$modelValue);
        }
    };
}]);

/***********************
 BUTTONS
 ***********************/

/**
 * @description Provides dropdown menu functionality in place of bootstrap js
 * @example
<li class="dropdown">
    <a class="dropdown-toggle">My Dropdown Menu</a>
    <ul class="dropdown-menu">
        <li ng-repeat="choice in dropChoices">
            <a ng-href="{{choice.href}}">{{choice.text}}</a>
        </li>
    </ul>
 </li>
 */
Application.Interface.directive('dropdownToggle', ['$document', '$location', function ($document, $location) {
    var openElement = null,
        closeMenu   = angular.noop;
    return {
        restrict: 'CA',
        link: function(scope, element, attrs) {
            scope.$watch('$location.path', function() { closeMenu(); });
            element.parent().bind('click', function() { closeMenu(); });
            element.bind('click', function (event) {

                var elementWasOpen = (element === openElement);

                event.preventDefault();
                event.stopPropagation();

                if (!!openElement) {
                    closeMenu();
                }

                if (!elementWasOpen) {
                    element.parent().addClass('open');
                    openElement = element;
                    closeMenu = function (event) {
                        if (event) {
                            event.preventDefault();
                            event.stopPropagation();
                        }
                        $document.unbind('click', closeMenu);
                        element.parent().removeClass('open');
                        closeMenu = angular.noop;
                        openElement = null;
                    };
                    $document.bind('click', closeMenu);
                }
            });
        }
    };
}]);

//from Angular UI Bootstrap
Application.Interface.directive('btnRadio', function () {
    var activeClass = 'active';
    var toggleEvent = 'click';

    return {
        require:'ngModel',
        link:function (scope, element, attrs, ngModelCtrl) {

            var value = scope.$eval(attrs.btnRadio);

            //model -> UI
            scope.$watch(function () {
                return ngModelCtrl.$modelValue;
            }, function (modelValue) {
                if (angular.equals(modelValue, value)){
                    element.addClass(activeClass);
                } else {
                    element.removeClass(activeClass);
                }
            });

            //ui->model
            element.bind(toggleEvent, function () {
                if (!element.hasClass(activeClass)) {
                    scope.$apply(function () {
                        ngModelCtrl.$setViewValue(value);
                    });
                }
            });
        }
    };
});

//from Angular UI Bootstrap
Application.Interface.directive('btnCheckbox', function() {

    var activeClass = 'active',
        toggleEvent = 'click';

    return {
        require:'ngModel',
        link:function (scope, element, attrs, ngModelCtrl) {

            var trueValue = scope.$eval(attrs.btnCheckboxTrue);
            var falseValue = scope.$eval(attrs.btnCheckboxFalse);

            trueValue = angular.isDefined(trueValue) ? trueValue : true;
            falseValue = angular.isDefined(falseValue) ? falseValue : false;

            //model -> UI
            scope.$watch(function () {
                return ngModelCtrl.$modelValue;
            }, function (modelValue) {
                if (angular.equals(modelValue, trueValue)) {
                    element.addClass(activeClass);
                } else {
                    element.removeClass(activeClass);
                }
            });

            //ui->model
            element.bind(toggleEvent, function () {
                scope.$apply(function () {
                    ngModelCtrl.$setViewValue(element.hasClass(activeClass) ? falseValue : trueValue);
                });
            });
        }
    };
});


/***********************
 KEYSTROKES
 ***********************/

/**
 * Bind one or more handlers to particular keys or their combination
 * @param hash {mixed} keyBindings Can be an object or string where keybinding expression of keys or keys combinations and AngularJS Exspressions are set. Object syntax: "{ keys1: expression1 [, keys2: expression2 [ , ... ]]}". String syntax: ""expression1 on keys1 [ and expression2 on keys2 [ and ... ]]"". Expression is an AngularJS Expression, and key(s) are dash-separated combinations of keys and modifiers (one or many, if any. Order does not matter). Supported modifiers are 'ctrl', 'shift', 'alt' and key can be used either via its keyCode (13 for Return) or name. Named keys are 'backspace', 'tab', 'enter', 'esc', 'space', 'pageup', 'pagedown', 'end', 'home', 'left', 'up', 'right', 'down', 'insert', 'delete', 'period', 'comma'.
 * @example <input ui-keypress="{enter:'x = 1', 'ctrl-shift-space':'foo()', 'shift-13':'bar()'}" /> <input ui-keypress="foo = 2 on ctrl-13 and bar('hello') on shift-esc" />
 * @note keypress for arrows and some other keys will not work -- use keydown or keyup in that case
 **/
Application.Interface.directive('uiKeydown', ['keypressHelper', function(keypressHelper){
    return {
        link: function (scope, elm, attrs) {
            keypressHelper('keydown', scope, elm, attrs);
        }
    };
}]);

Application.Interface.directive('uiKeypress', ['keypressHelper', function(keypressHelper){
    return {
        link: function (scope, elm, attrs) {
            keypressHelper('keypress', scope, elm, attrs);
        }
    };
}]);

Application.Interface.directive('uiKeyup', ['keypressHelper', function(keypressHelper){
    return {
        link: function (scope, elm, attrs) {
            keypressHelper('keyup', scope, elm, attrs);
        }
    };
}]);

//alternative to ng-pattern, which will allow all input but if model not acceptable then not set -- this will not allow invalid values to propagate to model, or be visible in the input field
Application.Interface.directive('restrictInput', ['$parse', function($parse) {
    return {
        restrict: 'A',
        require: 'ngModel',
        link: function(scope, iElement, iAttrs, ngModel) {
            var fn = function(input) {
                //in function so remains dynamic
                var regexp = scope.$eval(iAttrs.restrictInput);

                //testing
                //console.log("RESTRICT", ngModel.$viewValue, ngModel.$modelValue, input, ngModel);

                var transformedInput = input.replace(regexp, '');
                if (transformedInput != input) {
                    ngModel.$setViewValue(transformedInput);
                    ngModel.$render();
                    //testing
                    //console.log("RESTRICT 2", ngModel.$viewValue, ngModel.$modelValue, input, ngModel);


                }
                return transformedInput;
            };

            ngModel.$parsers.unshift( fn );
        }
    }
}]);

/***********************
 DIALOG / MODAL
 ***********************/

Application.Interface.directive('modal', ['$parse', '$dialog', function($parse, $dialog) {
    return {
        restrict: 'EA',
        terminal: true,
        link: function(scope, elm, attrs) {
            var opts = angular.extend({}, scope.$eval(attrs.options));
            var shownExpr = attrs.modal || attrs.show;
            var setClosed;

            // Create a dialog with the template as the contents of the directive
            // Add the current scope as the resolve in order to make the directive scope as a dialog controller scope
            opts = angular.extend(opts, {
                template: elm.html(),
                resolve: { $scope: function() { return scope; } }
            });
            var dialog = $dialog.dialog(opts);

            elm.remove();

            if (attrs.close) {
                setClosed = function() {
                    $parse(attrs.close)(scope);
                };
            } else {
                setClosed = function() {
                    if (angular.isFunction($parse(shownExpr).assign)) {
                        $parse(shownExpr).assign(scope, false);
                    }
                };
            }

            scope.$watch(shownExpr, function(isShown, oldShown) {
                if (isShown) {
                    dialog.open().then(function(){
                        setClosed();
                    });
                } else {
                    //Make sure it is not opened
                    if (dialog.isOpen()){
                        dialog.close();
                    }
                }
            });
        }
    };
}]);


/******************
 TOOLTIP + POPOVER
 ******************/

/**
 * The following features are still outstanding:
 * animation as a function
 * placement as a function
 * inside
 * support for more triggers than just mouse enter/leave
 * html tooltips
 * and selector delegation.
 */

/**
 * The $tooltip service creates tooltip- and popover-like directives as well as
 * houses global options for them.
 */
Application.Interface.provider( '$tooltip', function () {
    // The default options tooltip and popover.
    var defaultOptions = {
        placement: 'top',
        animation: true,
        popupDelay: 0
    };

    // Default hide triggers for each show trigger
    var triggerMap = {
        'mouseenter': 'mouseleave',
        'click': 'click',
        'focus': 'blur'
    };

    // The options specified to the provider globally.
    var globalOptions = {};

    /**
     * `options({})` allows global configuration of all tooltips in the
     * application.
     *
     *   var app = angular.module( 'App', ['ui.bootstrap.tooltip'], function( $tooltipProvider ) {
   *     // place tooltips left instead of top by default
   *     $tooltipProvider.options( { placement: 'left' } );
   *   });
     */
    this.options = function( value ) {
        angular.extend( globalOptions, value );
    };

    /**
     * This allows you to extend the set of trigger mappings available. E.g.:
     *
     *   $tooltipProvider.setTriggers( 'openTrigger': 'closeTrigger' );
     */
    this.setTriggers = function setTriggers ( triggers ) {
        angular.extend( triggerMap, triggers );
    };

    /**
     * This is a helper function for translating camel-case to snake-case.
     */
    function snake_case(name){
        var regexp = /[A-Z]/g;
        var separator = '-';
        return name.replace(regexp, function(letter, pos) {
            return (pos ? separator : '') + letter.toLowerCase();
        });
    }

    /**
     * Returns the actual instance of the $tooltip service.
     */
    this.$get = [ '$window', '$compile', '$timeout', '$parse', '$document', '$position', '$interpolate', function ( $window, $compile, $timeout, $parse, $document, $position, $interpolate ) {
        return function $tooltip ( type, prefix, defaultTriggerShow ) {
            var options = angular.extend( {}, defaultOptions, globalOptions );

            /**
             * Returns an object of show and hide triggers.
             *
             * If a trigger is supplied,
             * it is used to show the tooltip; otherwise, it will use the `trigger`
             * option passed to the `$tooltipProvider.options` method; else it will
             * default to the trigger supplied to this directive factory.
             *
             * The hide trigger is based on the show trigger. If the `trigger` option
             * was passed to the `$tooltipProvider.options` method, it will use the
             * mapped trigger from `triggerMap` or the passed trigger if the map is
             * undefined; otherwise, it uses the `triggerMap` value of the show
             * trigger; else it will just use the show trigger.
             */
            function setTriggers ( trigger ) {
                var show, hide;

                show = trigger || options.trigger || defaultTriggerShow;
                if ( angular.isDefined ( options.trigger ) ) {
                    hide = triggerMap[options.trigger] || show;
                } else {
                    hide = triggerMap[show] || show;
                }

                return {
                    show: show,
                    hide: hide
                };
            }

            var directiveName = snake_case( type );
            var triggers = setTriggers( undefined );

            var startSym = $interpolate.startSymbol();
            var endSym = $interpolate.endSymbol();
            var template =
                '<'+ directiveName +'-popup '+
                    'title="'+startSym+'tt_title'+endSym+'" '+
                    'content="'+startSym+'tt_content'+endSym+'" '+
                    'placement="'+startSym+'tt_placement'+endSym+'" '+
                    'animation="tt_animation()" '+
                    'is-open="tt_isOpen"'+
                    '>'+
                    '</'+ directiveName +'-popup>';

            return {
                restrict: 'EA',
                scope: true,
                link: function link ( scope, element, attrs ) {
                    var tooltip = $compile( template )( scope );
                    var transitionTimeout;
                    var popupTimeout;
                    var $body;
                    var appendToBody = angular.isDefined( options.appendToBody ) ? options.appendToBody : false;

                    // By default, the tooltip is not open.
                    // TODO add ability to start tooltip opened
                    scope.tt_isOpen = false;

                    function toggleTooltipBind () {
                        if ( ! scope.tt_isOpen ) {
                            showTooltipBind();
                        } else {
                            hideTooltipBind();
                        }
                    }

                    // Show the tooltip with delay if specified, otherwise show it immediately
                    function showTooltipBind() {
                        if ( scope.tt_popupDelay ) {
                            popupTimeout = $timeout( show, scope.tt_popupDelay );
                        } else {
                            scope.$apply( show );
                        }
                    }

                    function hideTooltipBind () {
                        scope.$apply(function () {
                            hide();
                        });
                    }

                    // Show the tooltip popup element.
                    function show() {
                        var position,
                            ttWidth,
                            ttHeight,
                            ttPosition;

                        // Don't show empty tooltips.
                        if ( ! scope.tt_content ) {
                            return;
                        }

                        // If there is a pending remove transition, we must cancel it, lest the
                        // tooltip be mysteriously removed.
                        if ( transitionTimeout ) {
                            $timeout.cancel( transitionTimeout );
                        }

                        // Set the initial positioning.
                        tooltip.css({ top: 0, left: 0, display: 'block' });

                        // Now we add it to the DOM because need some info about it. But it's not
                        // visible yet anyway.
                        if ( appendToBody ) {
                            $body = $body || $document.find( 'body' );
                            $body.append( tooltip );
                        } else {
                            element.after( tooltip );
                        }

                        // Get the position of the directive element.
                        position = appendToBody ? $position.offset( element ) : $position.position( element );

                        // Get the height and width of the tooltip so we can center it.
                        ttWidth = tooltip.prop( 'offsetWidth' );
                        ttHeight = tooltip.prop( 'offsetHeight' );

                        // Calculate the tooltip's top and left coordinates to center it with
                        // this directive.
                        switch ( scope.tt_placement ) {
                            case 'mouse':
                                var mousePos = $position.mouse();
                                ttPosition = {
                                    top: mousePos.y,
                                    left: mousePos.x
                                };
                                break;
                            case 'right':
                                ttPosition = {
                                    top: position.top + position.height / 2 - ttHeight / 2,
                                    left: position.left + position.width
                                };
                                break;
                            case 'bottom':
                                ttPosition = {
                                    top: position.top + position.height,
                                    left: position.left + position.width / 2 - ttWidth / 2
                                };
                                break;
                            case 'left':
                                ttPosition = {
                                    top: position.top + position.height / 2 - ttHeight / 2,
                                    left: position.left - ttWidth
                                };
                                break;
                            default:
                                ttPosition = {
                                    top: position.top - ttHeight,
                                    left: position.left + position.width / 2 - ttWidth / 2
                                };
                                break;
                        }

                        ttPosition.top += 'px';
                        ttPosition.left += 'px';

                        // Now set the calculated positioning.
                        tooltip.css( ttPosition );

                        //CUSTOM - avoid pointer events on tooltip - custom class
                        tooltip.addClass('noPointerEvents');

                        // And show the tooltip.
                        scope.tt_isOpen = true;
                    }

                    // Hide the tooltip popup element.
                    function hide() {
                        // First things first: we don't show it anymore.
                        scope.tt_isOpen = false;

                        //if tooltip is going to be shown after delay, we must cancel this
                        $timeout.cancel( popupTimeout );

                        // And now we remove it from the DOM. However, if we have animation, we
                        // need to wait for it to expire beforehand.
                        // FIXME: this is a placeholder for a port of the transitions library.
                        if ( angular.isDefined( scope.tt_animation ) && scope.tt_animation() ) {
                            transitionTimeout = $timeout( function () { tooltip.remove(); }, 500 );
                        } else {
                            tooltip.remove();
                        }
                    }

                    /**
                     * Observe the relevant attributes.
                     */
                    attrs.$observe( type, function ( val ) {
                        scope.tt_content = val;
                    });

                    attrs.$observe( prefix+'Title', function ( val ) {
                        scope.tt_title = val;
                    });

                    attrs.$observe( prefix+'Placement', function ( val ) {
                        scope.tt_placement = angular.isDefined( val ) ? val : options.placement;
                    });

                    attrs.$observe( prefix+'Animation', function ( val ) {
                        scope.tt_animation = angular.isDefined( val ) ? $parse( val ) : function(){ return options.animation; };
                    });

                    attrs.$observe( prefix+'PopupDelay', function ( val ) {
                        var delay = parseInt( val, 10 );
                        scope.tt_popupDelay = ! isNaN(delay) ? delay : options.popupDelay;
                    });

                    attrs.$observe( prefix+'Trigger', function ( val ) {
                        element.unbind( triggers.show );
                        element.unbind( triggers.hide );

                        triggers = setTriggers( val );

                        if ( triggers.show === triggers.hide ) {
                            element.bind( triggers.show, toggleTooltipBind );
                        } else {
                            element.bind( triggers.show, showTooltipBind );
                            element.bind( triggers.hide, hideTooltipBind );
                        }
                    });

                    attrs.$observe( prefix+'AppendToBody', function ( val ) {
                        appendToBody = angular.isDefined( val ) ? $parse( val )( scope ) : appendToBody;
                    });

                    // if a tooltip is attached to <body> we need to remove it on
                    // location change as its parent scope will probably not be destroyed
                    // by the change.
                    if ( appendToBody ) {
                        scope.$on('$locationChangeSuccess', function closeTooltipOnLocationChangeSuccess () {
                            if ( scope.tt_isOpen ) {
                                hide();
                            }
                        });
                    }

                    // Make sure tooltip is destroyed and removed.
                    scope.$on('$destroy', function onDestroyTooltip() {
                        if ( scope.tt_isOpen ) {
                            hide();
                        } else {
                            tooltip.remove();
                        }
                    });
                }
            };
        };
    }];
});

Application.Interface.directive( 'tooltipPopup', function () {
    return {
        restrict: 'E',
        replace: true,
        scope: { content: '@', placement: '@', animation: '&', isOpen: '&' },
        templateUrl: 'interface/templates/tooltip-popup.html'
    };
});

/**
 * @example
 * <span
 *      tooltip="{{dynamicTooltip}}"
 *      tooltip-placement="left"
 *      tooltip-animation="false"
 *      tooltip-popup-delay='1000'
 *      tooltip-trigger="focus"
 *      tooltip-html-unsafe="{{htmlTooltip}}" >
 * {{dynamicTooltipText}}
 * </span>
 */

Application.Interface.directive( 'tooltip', [ '$tooltip', function ( $tooltip ) {
    return $tooltip( 'tooltip', 'tooltip', 'mouseenter' );
}]);

Application.Interface.directive( 'tooltipHtmlUnsafePopup', function () {
    return {
        restrict: 'E',
        replace: true,
        scope: { content: '@', placement: '@', animation: '&', isOpen: '&' },
        templateUrl: 'interface/templates/tooltip-html-unsafe-popup.html'
    };
});

Application.Interface.directive( 'tooltipHtmlUnsafe', [ '$tooltip', function ( $tooltip ) {
    return $tooltip( 'tooltipHtmlUnsafe', 'tooltip', 'mouseenter' );
}]);


/**
 * The following features are still outstanding: popup delay, animation as a
 * function, placement as a function, inside, support for more triggers than
 * just mouse enter/leave, html popovers, and selector delegatation.
 */
Application.Interface.directive( 'popoverPopup', function () {
    return {
        restrict: 'EA',
        replace: true,
        scope: { title: '@', content: '@', placement: '@', animation: '&', isOpen: '&' },
        templateUrl: 'interface/templates/popover.html'
    };
});

/**
 * @example
 * <button class="btn"
 *  popover="{{dynamicPopover}}"
 *  popover-title="{{dynamicPopoverTitle}}"
 *  popover-placement="mouse"
 *  popover-trigger="mouseenter"
 *  popover-animation="true"
 *  popover-append-to-body="false" >
 */
Application.Interface.directive( 'popover', [ '$compile', '$timeout', '$parse', '$window', '$tooltip', function ( $compile, $timeout, $parse, $window, $tooltip ) {
        return $tooltip( 'popover', 'popover', 'click' );
    }]);

/***************
    DRAG N DROP - work in progress
 ***************/

/*
 * angular-dragon-drop v0.0.1
 * (c) 2013 Brian Ford http://briantford.com
 * License: MIT
 *
 * todo - clicking sometimes causes deletion of element
 * note - replaces ngRepeat in form draggable="item in collection"
 */


Application.Interface.directive('uiDrag', function ($document, $compile) {

    var dragValue,
        dragOrigin,
        floaty;

    var drag = function (ev) {
        var x = ev.clientX,
            y = ev.clientY;

        floaty.css('left', x + 10 + 'px');
        floaty.css('top', y + 10 + 'px');
    };

    // todo - use class unselectable (in base.css)
    var disableSelect = function () {
        angular.element(document.body).addClass('unselectable');
    };

    var enableSelect = function () {
        angular.element(document.body).removeClass('unselectable');
    };

    return {
        restrict: 'A',
        terminal: true,
        link: function (scope, elt, attr) {

            // get the `thing in things` expression
            console.log(attr);
            console.log(elt);
            var expression = attr.uiDrag;
            var match = expression.match(/^\s*(.+)\s+in\s+(.*?)\s*$/);
            if (!match) {
                throw Error("Expected ngRepeat in form of '_item_ in _collection_' but got '" +
                    expression + "'.");
            }
            var lhs = match[1];
            var rhs = match[2];

            //todo -this sucks - rewrite
            // pull out the template to re-use. Improvised ng-transclude.
            var template = elt.html();
            elt.html('');
            var child = angular.element('<div ng-repeat="' + lhs + ' in ' + rhs + '">' + template + '</div>');
            elt.append(child);

            $compile(child)(scope);

            var spawnFloaty = function () {
                scope.$apply(function () {
                    floaty = angular.element('<div style="position: fixed;">' + template + '</div>');
                    var floatyScope = scope.$new();
                    floatyScope[lhs] = dragValue;
                    $compile(floaty)(floatyScope);
                    angular.element(document.body).append(floaty);
                });

                $document.bind('mousemove', drag);
            };

            var killFloaty = function () {
                $document.unbind('mousemove', drag);
                if (floaty) {
                    floaty.remove();
                    floaty = null;
                }
            };

            elt.bind('mousedown', function (ev) {
                if (dragValue) {
                    return;
                }
                scope.$apply(function () {
                    var targetScope = angular.element(ev.target).scope();
                    var value = dragValue = targetScope[lhs];
                    //console.log(value);
                    var list = scope.$eval(rhs);
                    dragOrigin = list;
                    list.splice(list.indexOf(value), 1);
                });
                disableSelect();
                spawnFloaty();
                drag(ev);

            });

            // handle something being dropped here
            elt.bind('mouseup', function (ev) {
                if (dragValue) {
                    scope.$apply(function () {
                        var list = scope.$eval(rhs);
                        list.push(dragValue);
                        dragValue = dragOrigin = null;
                    });
                }
                enableSelect();
                killFloaty();
            });

            // else, the event bubbles up to document
            $document.bind('mouseup', function (ev) {
                if (dragValue) {
                    scope.$apply(function () {
                        dragOrigin.push(dragValue);
                        dragValue = dragOrigin = null;
                    });
                    enableSelect();
                    killFloaty();
                }
            });

        }
    };
});

//note - 'draggable' is an HTML5 term - will be true or false when used natively - should use another attr
//see http://jsfiddle.net/ADukg/2516/
Application.Interface.directive('uiDraggable', ['$compile', '$document', function($compile, $document) {

    var dragEl,
        dragOrigin,
        originalPositionCss;

    var drag = function (ev) {
        var x = ev.clientX,
            y = ev.clientY;

        dragEl.css('left', x + 10 + 'px');
        dragEl.css('top', y + 10 + 'px');
    };

    // todo - use class unselectable (in base.css)?
    var disableSelect = function () {
        angular.element(document.body).css({
            '-moz-user-select': '-moz-none',
            '-khtml-user-select': 'none',
            '-webkit-user-select': 'none',
            '-ms-user-select': 'none',
            'user-select': 'none'
        });
    };

    var enableSelect = function () {
        angular.element(document.body).css({
            '-moz-user-select': '',
            '-khtml-user-select': '',
            '-webkit-user-select': '',
            '-ms-user-select': '',
            'user-select': ''
        });
    };

    return {
        restrict: 'A',
        terminal: false,
        link: function (scope, element, attrs) {
            console.log(element);

            var targets = $('body').find(attrs.target);
            var template = element.html();
            //todo - use
            originalPositionCss = element.css('position');


            //todo - check for siblings after
            var saveLocation = function(element) {
                var loc = {};
                loc.parent = element.parent();
                return(loc);
            };

            var restoreLocation = function(loc) {
                if (angular.isUndefined(loc) || loc == null ||loc == {}) return;

                console.log('restoring');

                (loc.parent).prepend(dragEl);
            };


            dragOrigin = saveLocation(element);


            var spawnFloaty = function () {
                element.css({'position' : 'fixed'});
                $document.bind('mousemove', drag);
                console.log(drag);
            };

            var killFloaty = function () {
                $document.unbind('mousemove', drag);
                console.log(drag);
                element.css({position: 'static'});
                restoreLocation(dragOrigin);
            };

            //events

            element.bind('mousedown', function (ev) {
                if (dragEl) {
                    return;
                }

                scope.$apply(function() {
                    dragEl = element;
                    disableSelect();
                    spawnFloaty();
                    drag(ev);
                });
            });

            // handle something being dropped here
            // todo - rewrite for target
            element.bind('mouseup', function (ev) {
                console.log('el mouseup');
                if (dragEl) {
                    console.log('el mouseup - yup dragvalue');
                    scope.$apply(function() {
                        enableSelect();
                        killFloaty();
                        dragEl = dragOrigin = null;
                    });
                }

            });

            // else, the event bubbles up to document
            $document.bind('mouseup', function (ev) {
                console.log('ducment call');
                if (dragEl) {
                    console.log('mouseup doc action');
                    scope.$apply(function () {
                        enableSelect();
                        killFloaty();
                        dragEl = dragOrigin = null;
                    });

                }
            });
        }
    }
}]);


/*********************
    TYPEAHEAD
*********************/



/**
 * A helper service that can parse typeahead's syntax (string provided by users)
 * Extracted to a separate service for ease of unit testing
 */
Application.Interface.factory('typeaheadParser', ['$parse', function ($parse) {

        //                      00000111000000000000022200000000000000003333333333333330000000000044000
        var TYPEAHEAD_REGEXP = /^\s*(.*?)(?:\s+as\s+(.*?))?\s+for\s+(?:([\$\w][\$\w\d]*))\s+in\s+(.*)$/;

        return {
            parse:function (input) {

                var match = input.match(TYPEAHEAD_REGEXP), modelMapper, viewMapper, source;
                if (!match) {
                    throw new Error(
                        "Expected typeahead specification in form of '_modelValue_ (as _label_)? for _item_ in _collection_'" +
                            " but got '" + input + "'.");
                }

                return {
                    itemName:match[3],
                    source:$parse(match[4]),
                    viewMapper:$parse(match[2] || match[1]),
                    modelMapper:$parse(match[1])
                };
            }
        };
    }]);

Application.Interface.directive('typeahead', ['$compile', '$parse', '$q', '$timeout', '$document', '$position', 'typeaheadParser', function ($compile, $parse, $q, $timeout, $document, $position, typeaheadParser) {

        var HOT_KEYS = [9, 13, 27, 38, 40];

        return {
            require:'ngModel',
            link:function (originalScope, element, attrs, modelCtrl) {

                //SUPPORTED ATTRIBUTES (OPTIONS)

                //minimal no of characters that needs to be entered before typeahead kicks-in
                var minSearch = originalScope.$eval(attrs.typeaheadMinLength) || 1;

                //minimal wait time after last character typed before typehead kicks-in
                var waitTime = originalScope.$eval(attrs.typeaheadWaitMs) || 0;

                //should it restrict model values to the ones selected from the popup only?
                var isEditable = originalScope.$eval(attrs.typeaheadEditable) !== false;

                //binding to a variable that indicates if matches are being retrieved asynchronously
                var isLoadingSetter = $parse(attrs.typeaheadLoading).assign || angular.noop;

                //a callback executed when a match is selected
                var onSelectCallback = $parse(attrs.typeaheadOnSelect);

                var inputFormatter = attrs.typeaheadInputFormatter ? $parse(attrs.typeaheadInputFormatter) : undefined;

                //INTERNAL VARIABLES

                //model setter executed upon match selection
                var $setModelValue = $parse(attrs.ngModel).assign;

                //expressions used by typeahead
                var parserResult = typeaheadParser.parse(attrs.typeahead);


                //pop-up element used to display matches
                var popUpEl = angular.element('<typeahead-popup></typeahead-popup>');
                popUpEl.attr({
                    matches: 'matches',
                    active: 'activeIdx',
                    select: 'select(activeIdx)',
                    query: 'query',
                    position: 'position'
                });
                //custom item template
                if (angular.isDefined(attrs.typeaheadTemplateUrl)) {
                    popUpEl.attr('template-url', attrs.typeaheadTemplateUrl);
                }

                //create a child scope for the typeahead directive so we are not polluting original scope
                //with typeahead-specific data (matches, query etc.)
                var scope = originalScope.$new();
                originalScope.$on('$destroy', function(){
                    scope.$destroy();
                });

                var resetMatches = function() {
                    scope.matches = [];
                    scope.activeIdx = -1;
                };

                var getMatchesAsync = function(inputValue) {

                    var locals = {$viewValue: inputValue};
                    isLoadingSetter(originalScope, true);
                    $q.when(parserResult.source(scope, locals)).then(function(matches) {

                        //it might happen that several async queries were in progress if a user were typing fast
                        //but we are interested only in responses that correspond to the current view value
                        if (inputValue === modelCtrl.$viewValue) {
                            if (matches.length > 0) {

                                scope.activeIdx = 0;
                                scope.matches.length = 0;

                                //transform labels
                                for(var i=0; i<matches.length; i++) {
                                    locals[parserResult.itemName] = matches[i];
                                    scope.matches.push({
                                        label: parserResult.viewMapper(scope, locals),
                                        model: matches[i]
                                    });
                                }

                                scope.query = inputValue;
                                //position pop-up with matches - we need to re-calculate its position each time we are opening a window
                                //with matches as a pop-up might be absolute-positioned and position of an input might have changed on a page
                                //due to other elements being rendered
                                scope.position = $position.position(element);
                                scope.position.top = scope.position.top + element.prop('offsetHeight');

                            } else {
                                resetMatches();
                            }
                            isLoadingSetter(originalScope, false);
                        }
                    }, function(){
                        resetMatches();
                        isLoadingSetter(originalScope, false);
                    });
                };

                resetMatches();

                //we need to propagate user's query so we can higlight matches
                scope.query = undefined;

                //Declare the timeout promise var outside the function scope so that stacked calls can be cancelled later
                var timeoutPromise;

                //plug into $parsers pipeline to open a typeahead on view changes initiated from DOM
                //$parsers kick-in on all the changes coming from the view as well as manually triggered by $setViewValue
                modelCtrl.$parsers.push(function (inputValue) {

                    resetMatches();
                    if (inputValue && inputValue.length >= minSearch) {
                        if (waitTime > 0) {
                            if (timeoutPromise) {
                                $timeout.cancel(timeoutPromise);//cancel previous timeout
                            }
                            timeoutPromise = $timeout(function () {
                                getMatchesAsync(inputValue);
                            }, waitTime);
                        } else {
                            getMatchesAsync(inputValue);
                        }
                    }

                    return isEditable ? inputValue : undefined;
                });

                modelCtrl.$formatters.push(function (modelValue) {

                    var candidateViewValue, emptyViewValue;
                    var locals = {};

                    if (inputFormatter) {

                        locals['$model'] = modelValue;
                        return inputFormatter(originalScope, locals);

                    } else {
                        locals[parserResult.itemName] = modelValue;

                        //it might happen that we don't have enough info to properly render input value
                        //we need to check for this situation and simply return model value if we can't apply custom formatting
                        candidateViewValue = parserResult.viewMapper(originalScope, locals);
                        emptyViewValue = parserResult.viewMapper(originalScope, {});

                        return candidateViewValue!== emptyViewValue ? candidateViewValue : modelValue;
                    }
                });

                scope.select = function (activeIdx) {
                    //called from within the $digest() cycle
                    var locals = {};
                    var model, item;

                    locals[parserResult.itemName] = item = scope.matches[activeIdx].model;
                    model = parserResult.modelMapper(originalScope, locals);
                    $setModelValue(originalScope, model);

                    onSelectCallback(originalScope, {
                        $item: item,
                        $model: model,
                        $label: parserResult.viewMapper(originalScope, locals)
                    });

                    //return focus to the input element if a mach was selected via a mouse click event
                    resetMatches();
                    element[0].focus();
                };

                //bind keyboard events: arrows up(38) / down(40), enter(13) and tab(9), esc(27)
                element.bind('keydown', function (evt) {

                    //typeahead is open and an "interesting" key was pressed
                    if (scope.matches.length === 0 || HOT_KEYS.indexOf(evt.which) === -1) {
                        return;
                    }

                    evt.preventDefault();

                    if (evt.which === 40) {
                        scope.activeIdx = (scope.activeIdx + 1) % scope.matches.length;
                        scope.$digest();

                    } else if (evt.which === 38) {
                        scope.activeIdx = (scope.activeIdx ? scope.activeIdx : scope.matches.length) - 1;
                        scope.$digest();

                    } else if (evt.which === 13 || evt.which === 9) {
                        scope.$apply(function () {
                            scope.select(scope.activeIdx);
                        });

                    } else if (evt.which === 27) {
                        evt.stopPropagation();

                        resetMatches();
                        scope.$digest();
                    }
                });

                $document.bind('click', function(){
                    resetMatches();
                    scope.$digest();
                });

                element.after($compile(popUpEl)(scope));
            }
        };

    }]);

Application.Interface.directive('typeaheadPopup', function () {
        return {
            restrict:'E',
            scope:{
                matches:'=',
                query:'=',
                active:'=',
                position:'=',
                select:'&'
            },
            replace:true,
            templateUrl:'/interface/templates/typeahead-popup.html',
            link:function (scope, element, attrs) {

                scope.templateUrl = attrs.templateUrl;

                scope.isOpen = function () {
                    return scope.matches.length > 0;
                };

                scope.isActive = function (matchIdx) {
                    return scope.active == matchIdx;
                };

                scope.selectActive = function (matchIdx) {
                    scope.active = matchIdx;
                };

                scope.selectMatch = function (activeIdx) {
                    scope.select({activeIdx:activeIdx});
                };
            }
        };
    });

Application.Interface.directive('typeaheadMatch', ['$http', '$templateCache', '$compile', '$parse', function ($http, $templateCache, $compile, $parse) {
        return {
            restrict:'E',
            scope:{
                index:'=',
                match:'=',
                query:'='
            },
            link:function (scope, element, attrs) {
                var tplUrl = $parse(attrs.templateUrl)(scope.$parent) || '/interface/templates/typeahead-match.html';
                $http.get(tplUrl, {cache: $templateCache}).success(function(tplContent){
                    element.replaceWith($compile(tplContent.trim())(scope));
                });
            }
        };
    }]);

Application.Interface.filter('typeaheadHighlight', function() {

        function escapeRegexp(queryToEscape) {
            return queryToEscape.replace(/([.?*+^$[\]\\(){}|-])/g, "\\$1");
        }

        return function(matchItem, query) {
            return query ? matchItem.replace(new RegExp(escapeRegexp(query), 'gi'), '<strong>$&</strong>') : query;
        };
    });