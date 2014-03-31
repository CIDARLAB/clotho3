describe('clotho.core Collector', function() {
	var Collector,
		clothoLocalStorage,
		$window,
		prefix,
		mockStorage;

	beforeEach(function () {
		module('clotho.core');

		inject(function (_Collector_, _clothoLocalStorage_, _$window_) {
			Collector = _Collector_;
			clothoLocalStorage = _clothoLocalStorage_;
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
			console.log(key);
			delete mockStorage[key];
		});

		Object.defineProperty($window.localStorage, 'length', {
			get: function () {
				console.log('length got called: ' + Object.keys(this).length);
				return Object.keys(this).length - 2;
			}
		});

		spyOn($window.localStorage, 'key').and.callFake(function(index) {
			console.log('want key ' + index);
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


	it('should should return the whole collector object', function() {
		expect(angular.isObject(Collector.collector)).toBe(true);
	});

	describe('#storeModel', function() {

		it('should have a storeModel function', function() {
			expect(angular.isFunction(Collector.storeModel)).toBe(true);
		});

		it('should add an item when calling storeModel and retrieve it with retrieveModel', function() {
			var myObj = {hi: 'there'};
			Collector.storeModel('12849124219', myObj);

			expect(Collector.retrieveModel('12849124219')).toEqual(myObj);
		});

	});


	describe('#retrieveModel', function() {

		it('should have a retrieveModel function', function() {
			expect(angular.isFunction(Collector.retrieveModel)).toBe(true);
		});


		it('should prefix a key with clotho_', function() {
			var model = {myData : "is so hot"};
			Collector.storeModel('123', model);
			var direct = $window.localStorage.getItem('clotho_123');
			var indirect = Collector.retrieveModel('123');
			expect(angular.isObject(indirect)).toBe(true);
			expect(angular.isString(direct)).toBe(true);
			expect(JSON.stringify(indirect)).toEqual(direct);
		});

	});
});