'use strict';

// note - requires jQuery
// future - to make more universal, see:
    // https://raw.github.com/cowboy/jquery-outside-events/v1.1/jquery.ba-outside-events.js
    // angular way of creating ng-directives
Application.Interface.directive('clickOutside', ['$document', '$parse', function($document, $parse) {
    return function(scope, element, attr, ctrl) {
        var handler = function(event) {
            //todo - rewrite without jQuery
            if (element.has(event.target).length == 0)
                scope.$apply(
                    $parse(attr['clickOutside'])(scope, {$event: event})
                );
        };

        $document.bind('click', handler);
        scope.$on('$destroy', function() {
            $document.unbind('click', handler);
        });
    }
}]);

//todo - write - see angularStrap and AngularUI, bootstrapJS
//note - jQuery dependence
Application.Interface.directive('modal', ['Clotho', '$http', '$parse', '$templateCache', function(Clotho, $http, $parse, $templateCache) {
    return function(scope, element, attr, ctrl) {
        var background = $('<div class="modal-backdrop"></div>');

        //todo .... ummm?? async possible?
        var url = $parse(attr['modal'])(scope, {$event: event});
        scope.content = $http.get(url, {cache: $templateCache});

        element.bind('click', function(event) {
            scope.show(event);
        });

        scope.show = function(event) {
            background.appendTo('body').click(function(event) {
                scope.hide(event);
            });


        };

        scope.hide = function(event) {
            background.remove();
            event.preventDefault();
        }

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

//future - remove polyfill once available-- currently not present in angular (1.1.4)
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

//future - remove polyfill once available-- currently not present in angular (1.1.4)
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

//fixme - not working because ng-model required
Application.Interface.directive('btnSelectable', ['$parse', function($parse) {

    var activeClass = 'active',
        toggleEvent = 'click';

    return {
        //require:'ngModel',
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
}]);


Application.Interface.factory('keypressHelper', ['$parse', function keypress($parse){
    var keysByCode = {
        8: 'backspace',
        9: 'tab',
        13: 'enter',
        27: 'esc',
        32: 'space',
        33: 'pageup',
        34: 'pagedown',
        35: 'end',
        36: 'home',
        37: 'left',
        38: 'up',
        39: 'right',
        40: 'down',
        45: 'insert',
        46: 'delete'
    };

    var capitaliseFirstLetter = function (string) {
        return string.charAt(0).toUpperCase() + string.slice(1);
    };

    return function(mode, scope, elm, attrs) {
        var params, combinations = [];
        params = scope.$eval(attrs['ui'+capitaliseFirstLetter(mode)]);

        // Prepare combinations for simple checking
        angular.forEach(params, function (v, k) {
            var combination, expression;
            expression = $parse(v);

            angular.forEach(k.split(' '), function(variation) {
                combination = {
                    expression: expression,
                    keys: {}
                };
                angular.forEach(variation.split('-'), function (value) {
                    combination.keys[value] = true;
                });
                combinations.push(combination);
            });
        });

        // Check only matching of pressed keys one of the conditions
        elm.bind(mode, function (event) {
            // No need to do that inside the cycle
            var metaPressed = !!(event.metaKey && !event.ctrlKey);
            var altPressed = !!event.altKey;
            var ctrlPressed = !!event.ctrlKey;
            var shiftPressed = !!event.shiftKey;
            var keyCode = event.keyCode;

            // normalize keycodes
            if (mode === 'keypress' && !shiftPressed && keyCode >= 97 && keyCode <= 122) {
                keyCode = keyCode - 32;
            }

            // Iterate over prepared combinations
            angular.forEach(combinations, function (combination) {

                var mainKeyPressed = combination.keys[keysByCode[event.keyCode]] || combination.keys[event.keyCode.toString()];

                var metaRequired = !!combination.keys.meta;
                var altRequired = !!combination.keys.alt;
                var ctrlRequired = !!combination.keys.ctrl;
                var shiftRequired = !!combination.keys.shift;

                if (
                    mainKeyPressed &&
                        ( metaRequired === metaPressed ) &&
                        ( altRequired === altPressed ) &&
                        ( ctrlRequired === ctrlPressed ) &&
                        ( shiftRequired === shiftPressed )
                    ) {
                    // Run the function
                    scope.$apply(function () {
                        combination.expression(scope, { '$event': event });
                    });
                }
            });
        });
    };
}]);

/**
 * Bind one or more handlers to particular keys or their combination
 * @param hash {mixed} keyBindings Can be an object or string where keybinding expression of keys or keys combinations and AngularJS Exspressions are set. Object syntax: "{ keys1: expression1 [, keys2: expression2 [ , ... ]]}". String syntax: ""expression1 on keys1 [ and expression2 on keys2 [ and ... ]]"". Expression is an AngularJS Expression, and key(s) are dash-separated combinations of keys and modifiers (one or many, if any. Order does not matter). Supported modifiers are 'ctrl', 'shift', 'alt' and key can be used either via its keyCode (13 for Return) or name. Named keys are 'backspace', 'tab', 'enter', 'esc', 'space', 'pageup', 'pagedown', 'end', 'home', 'left', 'up', 'right', 'down', 'insert', 'delete'.
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