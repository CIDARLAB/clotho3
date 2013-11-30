//todo - needs major update if going to use
//todo - handle val() if present, default to text() --- has been done in another custom directive
//note - can just use ng-bind when not editing

angular.module('clotho.interface').directive('contenteditable', function($timeout) {

	return {
		require: '?ngModel',
		link: function(scope, element, attrs, ngModel) {

			//NGMODEL VERSION
			if(!ngModel) return; // do nothing if no ng-model

			// Specify how UI should be updated
			ngModel.$render = function() {
				//console.log(ngModel);
				var toShow = ngModel.$modelValue || '';
				if (angular.isObject(toShow)) {
					toShow = JSON.stringify(toShow)
				}
				element.html(toShow);
			};

			// Write data to the model
			function read() {
				var text = element.text();
				// When we clear the content editable the browser leaves a <br> behind
				// Unless no-strip-br attribute is provided then we strip this out
				if( !attrs.noStripBr && text == '<br>' ) {
					text = '';
				}
				if (text != ngModel.$modelValue) {
					//console.log('new text: ', text);
					ngModel.$setViewValue(text);
				}
				if (text === '') {
					// the cursor disappears if the contents is empty so refocus
					$timeout(function(){
						element.blur();
						element.focus();
					})
				}
			}

			// Listen for change events to enable binding
			element.on('input', function() {
				scope.$apply(read);
			});

			//read(); // initialize

		}
	};
});