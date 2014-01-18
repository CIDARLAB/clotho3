angular.module('clotho.interface').directive('jqDroppable', function() {
	return {
		restrict: 'A',
		require: ['?ngModel', '^?jqDraggable'],
		link: function(scope, element, attrs, ngModel) {

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
				drop: null,
				over: null,
				out: null,
				activate: null,
				inactivate: null
			};


			//model handling
			//todo - decide
			//reference : https://github.com/codef0rmer/angular-dragdrop
			if (ngModel) {

				ngModel.$render = function() {
					element.droppable("refresh")
				};

				callbacks.over = function(e, ui) {

				};

				callbacks.out = function(e, ui) {

				};

				callbacks.drop = function(e, ui) {

				};
			}


			//watchers
			//attr level listener for disabled state
			scope.$watch(function() {
					return scope.$eval(attrs.jqDroppableDisabled)
				},
				function(newval, oldval) {
					element.droppable({disabled: newval})
				});

			//custom listeners and jQuery UI Droppable options
			scope.$watch(attrs.jqDroppable, function (newval, oldval) {
				angular.forEach(newval, function (value, key) {
					if( callbacks[key] ){
						// wrap the callback
						value = combineCallbacks( callbacks[key], value );
					}
					element.droppable('option', key, value);
				})
			}, true);

			angular.forEach(callbacks, function(value, key ){
				opts[key] = combineCallbacks(value, opts[key]);
			});

			//finally, create
			element.droppable(opts);
		}
	}
});