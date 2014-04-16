describe('clotho.core clothoLocalStorage', function() {
	var clothoLocalStorage,
		$window,
		PubSub,
		prefix,
		mockStorage;

	beforeEach(function () {
		module('clotho.core');

		inject(function (_clothoLocalStorage_, _PubSub_, _$window_) {
			clothoLocalStorage = _clothoLocalStorage_;
			PubSub = _PubSub_;
			$window = _$window_;
		});

		//this mocks localStorage
		prefix = "clotho_";
		mockStorage = {
			"nonClothoItem" : "some Value for this item",
			"clotho_123" : JSON.stringify({"hey" : "there"}),
			"clotho_456" : JSON.stringify({"tester" : "object"})
		};

		spyOn($window.localStorage, 'getItem').and.callFake(function(key) {
			return mockStorage[key];
		});

		spyOn($window.localStorage, 'setItem').and.callFake(function(key, value) {
			return mockStorage[key] = value;
		});

		spyOn($window.localStorage, 'removeItem').and.callFake(function(key) {
			delete mockStorage[key];
		});

		Object.defineProperty($window.localStorage, 'length', {
			get: function () { 
				return Object.keys(this).length - 2;
			}
		});

		spyOn($window.localStorage, 'key').and.callFake(function(index) {
			var ordKeys = Object.keys(mockStorage).sort();
			return ordKeys[index];
		});

		spyOn($window.localStorage, 'clear').and.callFake(function() {
			for (var key in mockStorage) {
				if (mockStorage.hasOwnProperty(key)) {
					delete mockStorage[key];
				}
			}
		});
	});

	it('should mirror the localStorage interface', function () {
		expect(angular.isFunction(clothoLocalStorage.getItem)).toBe(true);
		expect(angular.isFunction(clothoLocalStorage.hasItem)).toBe(true);
		expect(angular.isFunction(clothoLocalStorage.removeItem)).toBe(true);
		expect(angular.isFunction(clothoLocalStorage.setItem)).toBe(true);
		expect(angular.isFunction(clothoLocalStorage.clear)).toBe(true);
	});

	it('should check if localStorage is supported', function () {
		expect(clothoLocalStorage.isSupported()).toBe(true);
	});

	it('should prefix with "clotho_"', function () {
		expect(clothoLocalStorage.getPrefix()).toEqual(prefix);
	});

	it('should not store empty keys', function () {
		expect(clothoLocalStorage.setItem('', {some : "value"})).toBe(false);
	});

	it('should not store empty values', function () {
		expect(clothoLocalStorage.setItem('13833', '')).toBe(false);
	});

	it('should store strings as strings', function () {
		var exampleKey = '789';
		var exampleValue = 'here is a value';
		expect(clothoLocalStorage.setItem(exampleKey, exampleValue)).toBe(true);
		expect(mockStorage[prefix + exampleKey]).toEqual(exampleValue);
	});

	it('should store objects as strings', function () {
		var exampleKey = '987';
		var exampleValue = {some: "value"};
		expect(clothoLocalStorage.setItem(exampleKey, exampleValue)).toBe(true);
		expect(mockStorage[prefix + exampleKey]).toEqual(JSON.stringify(exampleValue));
	});

	it('should return objects', function () {
		expect(angular.isObject(clothoLocalStorage.getItem('123'))).toBe(true);
	});

	it('should return null if object not found', function () {
		expect(clothoLocalStorage.getItem('notPresent')).toBeNull();
	});

	it('should return default value when provided if object not found', function () {
		expect(clothoLocalStorage.getItem('notPresent', 'default')).toEqual('default');
	});

	it('should be able to check if an key exists', function () {
		expect(clothoLocalStorage.hasItem('123')).toBe(true);
		expect(clothoLocalStorage.hasItem('321')).toBe(false);
	});
	
	it('should not clear items without the prefix', function () {
		clothoLocalStorage.clear();
		expect(angular.isDefined(mockStorage.clotho_123)).toBe(false);
		expect(clothoLocalStorage.hasItem('123')).toBe(false);
		expect(angular.isDefined(mockStorage.nonClothoItem)).toBe(true);
	});

	//todo - test storage events

});