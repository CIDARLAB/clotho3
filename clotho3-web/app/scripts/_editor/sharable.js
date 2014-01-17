angular.module('clotho.editor').controller('Editor_SharableCtrl', function($scope, $compile, Clotho) {

	$scope.$watch('sharable', function(newval) {
		console.log(newval);
		Clotho.get(newval.schema).then(function(result) {
			$scope.schema = result;
			$scope.schema_custom = result.custom;

			//todo - this is kinda a hack... should get element more directly
			var insert = angular.element(document).find('insert-fields').html(generateDynamicForm());
			$compile(insert.contents())($scope);
		});
	});


	/**
	 * @description Generates HTML for a form provided a properly formed schema. Implicit parameters are scope, with $scope.schema defined
	 * @returns {string}
	 */
	function generateDynamicForm () {
		var fulltext = "";

		//todo - need BSON types from server (GH #99)
		/*
		 string -> test vs. textarea
		 when use select?
		 what to do for objects?
		 id not editable
		 */
		var typeMap = {
			'string' : 'text',
			'boolean' : 'checkbox'
		};

		angular.forEach($scope.schema.fields, function(field) {


			var type = field.type || 'text';
			if (type == '?') field.type == 'text';
			var required = field.required ? "required" : "";

			var htmlText_pre = '<form-field name="' + field.name + '">';
			var htmlText_post = '</form-field>';
			var inputText;

			switch (type) {
				case "textarea": {
					inputText = '<textarea rows="2" ' + required + ' ng-model="sharable.'+field.name+'"></textarea>';
					break;
				}
				case "select": {
					var optionsText = "";
					//todo - use ng-options + attach array to scope
					angular.forEach(field.options, function(value, key) {
						optionsText = optionsText + '<option value="'+value+'">'+ value + '</option>';
					});

					inputText = '<select ' + required + ' ng-model="sharable.'+field.name+'">' + optionsText + '</select>';
					break;
				}
				case "sharable": {
				}
				//todo - add filedrop support, and radio. checkbox works.
				default: {
					inputText = '<input type="' + type + '" ' + required + 'ng-model="sharable.'+field.name+'" >';
					break;
				}

			}

			fulltext += htmlText_pre + inputText + htmlText_post;
		});
		return fulltext;
	}


});