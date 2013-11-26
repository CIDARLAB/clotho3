angular.module('clotho.interface').directive('jqDraggable', function() {
	return {
		restrict : 'A',
		require: ['?ngModel', '^?jqDroppable'],
		compile: function compile(tElement, tAttrs, transclude) {
			return {
				pre: function preLink(scope, element, attrs) {
					scope.$watch(function() {
						return attrs.jqDraggableFloat
					}, function (newval, oldval) {
						if (!!newval)
							element.css({'position' : 'absolute'});
					});

				},
				post: function (scope, element, attrs, ngModel) {

					//set up
					function combineCallbacks(first,second){
						if( second && (typeof second === "function") ){
							return function(e,ui){
								first(e,ui);
								second(e,ui);
							};
						}
						return first;
					}

					var opts = {};
					var callbacks = {
						create: null,
						drag: null,
						start: null,
						stop: null
					};


					//model handling
					//todo - decide
					//reference : https://github.com/codef0rmer/angular-dragdrop
					if (ngModel) {
						console.log(ngModel);

						ngModel.$render = function() {
							element.draggable("refresh")
						};

						callbacks.start = function(e, ui) {

						};

						callbacks.stop = function(e, ui) {

						};

						callbacks.drag = function(e, ui) {

						};
					}


					//watchers
					//attr level listener for disabled state
					scope.$watch(function() {
							return scope.$eval(attrs.jqDraggableDisabled)
						},
						function(newval, oldval) {
							console.log(!!newval);
							element.draggable({disabled: !!newval})
						});

					//custom listeners and jQuery UI Draggable options
					scope.$watch(attrs.jqDraggable, function (newval, oldval) {
						angular.forEach(newval, function (value, key) {
							if( callbacks[key] ){
								// wrap the callback
								value = combineCallbacks( callbacks[key], value );
							}
							element.draggable('option', key, value);
						})
					}, true);

					angular.forEach(callbacks, function(value, key ){
						opts[key] = combineCallbacks(value, opts[key]);
					});

					//create
					element.draggable(opts);
				}
			}
		}
	}
});