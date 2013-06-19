'use strict';

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
Application.Interface.directive('clickOutside', ['$document', '$parse', function($document, $parse) {
    return function(scope, element, attr, ctrl) {
        var handler = function(event) {
            //todo - rewrite without jQuery.has()
            if (element.has(event.target).length == 0)
                scope.$safeApply(
                    $parse(attr['clickOutside'])(scope, {$event: event})
                );
        };

        $document.bind('click', handler);
        //todo - make sure this works
        scope.$on('$destroy', function() {
            $document.unbind('click', handler);
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
BUTTONS
 ***********************/

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